function coherency_output = perform_coherency(LFP,coh,my_experiment_variables)
% Function to estimate power coherency between two regions within the same
% phase
% Inputs:
% LFP: a structure containing the LFP data with the following
% fields:
% - fs: sampling frequency in Hz
% - type: type of the subject according to the subject file name (exp,ctrl,...)
% - LFP_data.(phase_name).(region_name): is a matrix of n_timepoints x n_channels 
%   that contains the LFP data of the region (i.e. PFC) of the particular
%   phase (i.e. baseline) 
% coh: a structure containing the parameters for the analyses, has the
% following fields:
%       - resolution: a positive number to determine the required frequency
%       resolution (0.1 or 1 for example)
%       - limits: an array of lenght two of positive numbers ([f1,f2], f2>f1) for frequency
%       limits to calculate coherency within that range (i.e. [0,100])
%       - pooling: 0 for no pooling, nan for automatic detection,other numbers for the minutes you want to pool
%       - window:  integer value used for windowing the data into that number of points to estimate the coherency 
%       - noverlap: integer value used for amount of overlapping datapoints between windows (must be smaller than window)
%       - ch_lvl: 0 for channel average, 1 for channel level analyses
% my_experiment_variable: the structure defined in the ephys_analyses
% routine containing parameters of the experiment
% Outputs:
% coherency_output : coherency_output.(phase_name).(str_region) is a structure that
% contains the coherency results for the phase = phase_name and region =
% region1_region2 (based on the my_experiment_variables, for example: PFC_rs) and contains the
% following fields:
%   - coherency: the calculated coherency mscohere function. the size
%   will be [f1:freq_resoloution:f2] x (n_trials), n_trials is
%   determined based on pooling and wheter channel level analyses was 
%   requested or not. for example for a data of
%   10 minutes length and 32 channels, if the pooling is 1 and the ch_lvl
%   is 1, it will be 10*32, take note that first there will be the pooled versions 
%   of the same channel and then the next channel. if the ch_lvl = 0, it
%   will first average data over channels, then calculates the coherency
%   - freq: an array of [f1:freq_resoloution:f2] x 1 which
%   contains the frequencies at which the coherency was calculated
% --------------------------------------------------------------------
% the magnitude square coherency is calculated by dividing the signal into
% number of points as requested coh.window with the overlapping number of
% points requested in coh.noverlap, after calculating the coherency in each
% segment, the final coherency is calculated by taking average over the the
% coherency in each segment (check mscohere.m function)

%% initialize 
regions = my_experiment_variables.regions;
n_regions = length(regions);
phases = my_experiment_variables.phase;
n_phases = length(phases);
fs = LFP.fs;
freq_res = coh.resolution;
win = coh.window;
noverlap = coh.noverlap;
freq_start = coh.limits(1);
freq_end = coh.limits(2);
freq = freq_start:freq_res:freq_end;
confidence = coh.confidence;

coherency_output = struct;
%% calculate
for phase=1:n_phases % loop over phase 
    for reg1=1:n_regions % loop over region 1
        data_1 = LFP.(phases{phase}).(regions{reg1});
        for reg2=reg1+1:n_regions % loop over region 2
            data_2 = LFP.(phases{phase}).(regions{reg2});
            if coh.ch_lvl == 0 % average over cahnnels if channel wise is not requested
                data_avg_1 = mean(data_1,2);
                data_avg_2 = mean(data_2,2);
            else
                data_avg_1 = data_1;
                data_avg_2 = data_2;
            end
            
            if coh.pooling == 0 % check pooling
                len_pool = [];
            elseif isnan(coh.pooling)
                len_pool_1 = automatic_pooling(data_avg_1,fs);
                len_pool_2 = automatic_pooling(data_avg_2,fs);
                len_pool = max(len_pool_1,len_pool_2); % get the maximum number for pooling
            else
                len_pool = coh.pooling*60*fs;
                if rem(size(data_avg_1,1),len_pool) ~= 0
                    error('the requested pooling length and data length are not compatible');
                end
            end
            
            n_ch1 = size(data_avg_1,2);
            n_ch2 = size(data_avg_2,2);
            count = 1; % counter for all possible combinations
            C = [];
            C_th = [];
            for ch1=1:n_ch1 % loop over channels in the first region
                data_ch1 = data_avg_1(:,ch1);
                for ch2=1:n_ch2
                    data_ch2 = data_avg_2(:,ch2); % loop over channels in the second region
                    if isempty(len_pool)
                        [C(:,count),f] = mscohere(data_ch1,data_ch2,win,noverlap,freq,fs);
                        count = count + 1;
                    else
                        data_pool_1 = reshape(data_ch1,len_pool,[]); % reshape data based on pooling
                        data_pool_2 = reshape(data_ch2,len_pool,[]);
                        n_trials = size(data_pool_1,2);
                        C_trial = [];
                        C_th_tr = [];
                        for tr=1:n_trials % loop over trials
                            data_trial_1 = data_pool_1(:,tr);
                            data_trial_2 = data_pool_2(:,tr);
                            [C_trial(:,tr),f_trial(:,tr)] = mscohere(data_trial_1,data_trial_2,win,noverlap,freq,fs);
                            if confidence
                                [~,C_th_tr(:,tr)] = cohere_bootstrap_signif_level(data_trial_1,data_trial_2,0.95,100,win,noverlap,freq,fs);
                            end
                            %count = count + 1;
                        end
                        if confidence
                            C_th(:,count:count+n_trials-1) = C_th_tr;
                        end
                        C(:,count:count+n_trials-1) = C_trial;
                        count = size(C,2) + 1; % determine the end of C
                        f = f_trial(:,end);
                    end
                end
            end
            region_1 = regions{reg1};
            region_2 = regions{reg2};
            str_region = strcat(region_1,'_',region_2); % make the string as region1_region2 (i.e. PFC_rs)
            coherency_output.(phases{phase}).(str_region).coherency = C;
            coherency_output.(phases{phase}).(str_region).freq = f;  
            if confidence
                coherency_output.(phases{phase}).(str_region).coh_th = C_th; 
            end
        end
    end  
end
disp('Coherency analysis terminated ')

end

