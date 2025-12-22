function GC_output = perform_GC(LFP,cause,my_experiment_variables)
% Function to estimate conditional granger causality
%   Detailed explanation goes here
%Locally Weighted linear regression performed on time series downsampled
%data

GC_output = struct;

regions = my_experiment_variables.regions;
n_regions = length(regions);
phases = my_experiment_variables.phase;
n_phases = length(phases);
fs = LFP.fs;
%Down samplign factor govern frequency resolution AND frequency limits
down_sampling_factor = cause.down_sampling_factor;
len_trial = cause.len_trial*60*fs;
f_res = cause.resolution;
order_to_consider = cause.max_order;
icregmode = cause.icregmode;
regmode = cause.regmode;
fs_downsampled = fs/down_sampling_factor;
freq_to_estimate = 0:f_res:fs_downsampled/2;
alpha     = cause.alpha;
mhtc = cause.mhtc;
n_shuffles = cause.n_shuffles;

if cause.pooling == 0
    len_pooling = [];
else
    len_pooling = cause.pooling*60*fs;
end

for phase=1:n_phases
    data = [];
    %Average data over channels and pool
    for reg=1:n_regions
        if cause.ch_lvl 
            data.(regions{reg}) = LFP.(phases{phase}).(regions{reg});
        else
             data.(regions{reg}) = mean(LFP.(phases{phase}).(regions{reg}),2);
        end
    end
    n_ch1 = size(data.(regions{1}),2);
    n_ch2 = size(data.(regions{2}),2);
    count = 1;
    for ch1=1:n_ch1
        for ch2 = 1:n_ch2
            ch = [ch1,ch2];
            data_pool = [];
            for reg=1:n_regions
                if isempty(len_pooling)
                    data_pool(reg,:) = data.(regions{reg})(:,ch(reg));
                else
                    if rem(size(data_avg,1),len_pooling) ~= 0
                        error('the requested pooling length and data length are not compatible');
                    end
                    data_pool(reg,:,:) = reshape(data_avg,len_pooling,[]);
                end
            end
            n_pools = size(data_pool,3);
            for pool=1:n_pools
                tmp_data = data_pool(:,:,pool);
                if rem(size(tmp_data,2),len_trial) ~= 0
                    error('the requested trial length and data length are not compatible');
                end
                data_trial = [];
                for reg=1:n_regions
                    tmp_data_trial = tmp_data(reg,:);
                    data_trial(reg,:,:) = reshape(tmp_data_trial,len_trial,[]);
                end
                detreanded_data = mvdetrend(data_trial); %remove continuous linear trend
                swaped_data = permute(detreanded_data,[2 1 3]); %data in matrix time x regions x trials
                X_downsampled = movmeandecimatrix(swaped_data,down_sampling_factor);
                downsampled_data = permute(X_downsampled,[2 1 3]);
                [nvar,num_downsampled_data ,n_trials] = size(downsampled_data);
                [~,~,moAIC,moBIC] = tsdata_to_infocrit(downsampled_data,order_to_consider,icregmode);
                mo = min([moAIC,moBIC]);        %Pick lowest model order
                [BMat,EMat,Eseries] = tsdata_to_var(downsampled_data,mo,regmode);
                assert(~isbad(BMat),'VAR estimation failed');
                [~,RSQADJ] = rsquared(downsampled_data,Eseries); %Compute adjusted R2 of the  best model
                GC_output.(phases{phase}).Rsquared(:,count) = RSQADJ;
                [GammaMat,GammaInfo] = var_to_autocov(BMat,EMat); %Get autocovariance lags sequence (cov(X_t,X_t-k)
                while GammaInfo.rho>0.99
                    %Reduce model order till spectral radius is less than 0.99
                    %(why?)
                    mo = mo-1;
                    [BMat,EMat,Eseries] = tsdata_to_var(downsampled_data,mo,regmode);
                    assert(~isbad(BMat),'VAR estimation failed');
                    [~,RSQADJ] = rsquared(downsampled_data,Eseries);
                    GC_output.(phases{phase}).Rsquared(:,count) = RSQADJ;
                    [GammaMat,GammaInfo] = var_to_autocov(BMat,EMat);
                end
                granger_time = autocov_to_pwcgc(GammaMat); 
                % compute time-domain pairwise-conditional multivariate GC.
                %Note that if we have two regions, this is the same as computing
                %the simple GC  
                pval = mvgc_pval(granger_time,mo,num_downsampled_data,n_trials,1,1,nvar-2,[]); % take careful note of arguments!
                sig  = significance(pval,alpha,mhtc);
                for source=1:n_regions
                    for target=source+1:n_regions
                        %Spectral domani GC
                        f_1 = autocov_to_smvgc(GammaMat,source,target); % from target to source 
                        f_2 = autocov_to_smvgc(GammaMat,target,source); % from source to target
                        frequency_estimated = linspace(0,fs_downsampled/2,length(f_1));
                        str_1 = strcat(regions{target},'_to_',regions{source});
                        str_2 = strcat(regions{source},'_to_',regions{target});
                        GC_output.(phases{phase}).(str_1).infreq(:,count) = abs(interp1(frequency_estimated,f_1,freq_to_estimate,'spline'));
                        GC_output.(phases{phase}).(str_2).infreq(:,count) = abs(interp1(frequency_estimated,f_2,freq_to_estimate,'spline'));
                        GC_output.(phases{phase}).(str_1).intime(:,count) = granger_time(source,target);
                        GC_output.(phases{phase}).(str_2).intime(:,count) = granger_time(target,source);
                        GC_output.(phases{phase}).(str_1).issig(:,count) = sig(source,target);   
                        GC_output.(phases{phase}).(str_2).issig(:,count) = sig(target,source);
                        GC_output.(phases{phase}).(str_1).freq = freq_to_estimate;
                        GC_output.(phases{phase}).(str_2).freq = freq_to_estimate;
                        f_1_sh_all = [];
                        f_2_sh_all = [];
                        %compute shuffled statistics to assess significance
                        for n_sh=1:n_shuffles
                            downsampled_data_sh = downsampled_data(:,randperm(num_downsampled_data),:);
                            [BMat_sh,EMat_sh] = tsdata_to_var(downsampled_data_sh,mo,regmode);
                            assert(~isbad(BMat_sh),'VAR estimation failed');
                            [GammaMat_sh] = var_to_autocov(BMat_sh,EMat_sh);
                            f_1_sh = autocov_to_smvgc(GammaMat_sh,source,target); % from target to source
                            f_2_sh = autocov_to_smvgc(GammaMat_sh,target,source); % from source to target
                            frequency_estimated = linspace(0,fs_downsampled/2,length(f_1_sh));
                            f_1_sh_all(n_sh,:) = abs(interp1(frequency_estimated,f_1_sh,freq_to_estimate,'spline'));
                            f_2_sh_all(n_sh,:) = abs(interp1(frequency_estimated,f_2_sh,freq_to_estimate,'spline'));
                        end
                        GC_output.(phases{phase}).(str_1).infreq_sh(:,count) = prctile(f_1_sh_all,95);
                        GC_output.(phases{phase}).(str_2).infreq_sh(:,count) = prctile(f_2_sh_all,95);
                    end
                    
                end
                count = count+1;                
            end
        end
    end
end

end


