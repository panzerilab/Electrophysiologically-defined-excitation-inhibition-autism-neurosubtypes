function BLP_output = perform_BLP(LFP,BLP,my_experiment_variables)
% Function to calculate Band limited power correlation between regions
% within the same phase
% Inputs:
% LFP: a structure containing the LFP data with the following
% fields:
% - fs: sampling frequency in Hz
% - type: type of the subject according to the subject file name (exp,ctrl,...)
% - LFP_data.(phase_name).(region_name): is a matrix of n_timepoints x n_channels 
%   that contains the LFP data of the region (i.e. PFC) of the particular
%   phase (i.e. baseline) 
% BLP: a structure containing the parameters for the analyses, has the
% following fields:
%   - filterorder: integer determining the butterwurth filter order to use
%   - frequency: an array of lenght two of positive numbers ([f1,f2], f2>f1>=0) for the bandpass filter (if u want infraslow, use 0 instead 0.1 as the low frequency)
%   - pooling: 0 for no pooling, nan for automatic detection,other numbers for the minutes you want to pool
% my_experiment_variable: the structure defined in the ephys_analyses
% routine containing parameters of the experiment
% Outputs:
% BLP_output : BLP_output.(phase_name).(str_region) is a structure that
% contains the BLP results for the phase = phase_name and region =
% region1_region2 (based on the my_experiment_variables, for example: PFC_rs) and contains the
% following fields:
%   -BLP: the calculated correlation value between bandlimited power of LFP in region 1 and LFP in region 2. the size
%         will be 1 x n_trials. the n_trials depends on the pooling. for example for a data of
%   10 minutes length, if the pooling is 1 minute the n_trials will be 10
% --------------------------------------------------------------------
% at each phase and each region, the routine bandpass filters the signal using a
% non-causal user specified (filter order and bandwidth) butterwurth
% filter, then calculates the envelope of the signal as the envelope of 
% the Hilbert transformed signal. Then for each channel in region 1 and
% each channel in region2 it calculates the correlation between the envelopes 
% then averages the extracted correlations over channels. if pooling is non-zero it 
% does the same over each segmented part 

%% initialize
regions = my_experiment_variables.regions;
n_regions = length(regions);
phases = my_experiment_variables.phase;
n_phases = length(phases);
fs = LFP.fs;
freq_band = BLP.frequency;
order = BLP.filterorder;
[b,a] = butter(order,freq_band./(fs/2),'bandpass'); % create the filter

freq_res = BLP.resolution;
win = BLP.window;
noverlap = BLP.noverlap;
freq_start = BLP.limits(1);
freq_end = BLP.limits(2);
freq = freq_start:freq_res:freq_end; %Resolution of coherence of BLP envelope

% %% calculate
% for phase=1:n_phases % loop over phase
%     for reg1=1:n_regions % loop over region 1
%         data_1 = LFP.(phases{phase}).(regions{reg1});
%         data_1_LP = filtfilt(b,a,data_1); % filter the signal to the desired band
%         env_sig_1 = envelope(data_1_LP); % extract the envelope
%         for reg2=reg1+1:n_regions % loop over region 2
%             data_2 = LFP.(phases{phase}).(regions{reg2});
%             data_2_LP = filtfilt(b,a,data_2); % filter the signal
%             env_sig_2 = envelope(data_2_LP); % exctract the envelope
%             n_ch1 = size(env_sig_1,2);
%             n_ch2 = size(env_sig_2,2);
%             count = 1;
%             BLP_corr = [];
%             for ch1=1:n_ch1 % loop over channels in region 1
%                 data_ch1 = env_sig_1(:,ch1);
%                 for ch2=1:n_ch2 % loop over channels in region 2
%                     data_ch2 = env_sig_2(:,ch2);
%                     if BLP.pooling == 0 % check for pooling
%                         data_to_corr_1 = data_ch1;
%                         data_to_corr_2 = data_ch2;
%                     elseif isnan(BLP.pooling)
%                         len_pool_1 = automatic_pooling(mean(data_1,2),fs);
%                         len_pool_2 = automatic_pooling(mean(data_2,2),fs);
%                         len_pool = max(len_pool_1,len_pool_2);
%                         data_to_corr_1 = reshape(data_ch1,len_pool,[]);
%                         data_to_corr_2 = reshape(data_ch2,len_pool,[]);
%                     else
%                         len_pool = BLP.pooling*60*fs;
%                         if rem(size(data_ch1,1),len_pool) == 0
%                             data_to_corr_1 = reshape(data_ch1,len_pool,[]);
%                             data_to_corr_2 = reshape(data_ch2,len_pool,[]);
%                         else
%                             error('the requested pooling length and data length are not compatible');
%                         end
%                     end
%                     n_trials = size(data_to_corr_1,2);
%                     for tr=1:n_trials % loop over trails
%                         data_pooled_1 = data_to_corr_1(:,tr);
%                         data_pooled_2 = data_to_corr_2(:,tr);
%                         %BLP_corr(count,tr) = corr(data_pooled_1,data_pooled_2); % calculate the correlation between signals
%                         [C,F] = mscohere(data_pooled_1,data_pooled_2,win,noverlap,freq,fs);
%                         BLP_corr(count,tr) = mean(C);
%                     end
%                     count = count + 1;
%                 end
%             end
%             BLP_final = mean(BLP_corr,1); % average over channels
%             region_1 = regions{reg1};
%             region_2 = regions{reg2};
%             str_region = strcat(region_1,'_',region_2); % make the string as region1_region2 (i.e. PFC_rs)
%             BLP_output.(phases{phase}).(str_region).BLP = BLP_final;
%         end
%     end
% end
[b,a] = butter(order,freq_band./(fs/2),'bandpass');
    function [envl,fs_envelope] = compute_envelope(data)
        %Compute envelope as a sum over time points of the spectrogram
        [S,f_spectro,t_spec] = pspectrum(data, fs, 'spectrogram','FrequencyLimits',[0 80],...
                    'timeResolution',2,'OverlapPercent',75,'Leakage',0.85);
        S_gamma = log10(S(f_spectro>freq_band(1) & f_spectro<=freq_band(2),:));
        envl = (sum(S_gamma)) - mean((sum(S_gamma)));
        fs_envelope = 1./(t_spec(2)-t_spec(1));
    end
    function [envl,fs_envelope] = compute_envelope2(data)
        %Compute envelope as a sum over time points of the spectrogram
        data_filt = filtfilt(b,a,data); % filter the signal to the desired band
        envl= envelope(data_filt); % extract the envelope
        %fs_envelope = 1./(t_spec(2)-t_spec(1));
    end

%% calculate as integrated power of the spectrogram
 for phase=1:n_phases % loop over phase
     for reg1=1:n_regions % loop over region 1
         data_1 = LFP.(phases{phase}).(regions{reg1});
         data_1_m = mean(data_1,2);
         for reg2=reg1+1:n_regions % loop over region 2
            data_2 = LFP.(phases{phase}).(regions{reg2});
            data_2_m = mean(data_2,2);
            if BLP.pooling == 0 % check for pooling
                data_1_reshaped = data_1_m;
                data_2_reshaped = data_2_m;
            elseif isnan(BLP.pooling)
                len_pool_1 = automatic_pooling(mean(data_1,2),fs);
                len_pool_2 = automatic_pooling(mean(data_2,2),fs);
                len_pool = max(len_pool_1,len_pool_2);
                data_1_reshaped = reshape(data_1_m,len_pool,[]);
                data_2_reshaped = reshape(data_2_m,len_pool,[]);
            else
                len_pool = BLP.pooling*60*fs;
                if rem(size(data_1_m,1),len_pool) == 0
                    data_1_reshaped = reshape(data_1_m,len_pool,[]);
                    data_2_reshaped = reshape(data_2_m,len_pool,[]);
                else
                    error('the requested pooling length and data length are not compatible');
                end
            end
            BLP_corr = [];
            n_trials = size(data_1_reshaped,2);
            for tr=1:n_trials % loop over trails
                [env_sig_1,fs_envelope] = compute_envelope(data_1_reshaped(:,tr));
                [env_sig_2,fs_envelope] = compute_envelope(data_2_reshaped(:,tr));
                [C(:,tr),fc_envelope] = mscohere(env_sig_1,env_sig_2,[],[],freq,fs_envelope);
                %BLP_corr(tr) = mean(C);
            end
            BLP_final = C; % average over channels
            region_1 = regions{reg1};
            region_2 = regions{reg2};
            str_region = strcat(region_1,'_',region_2); % make the string as region1_region2 (i.e. PFC_rs)
            BLP_output.(phases{phase}).(str_region).fc_envelope = fc_envelope;
            BLP_output.(phases{phase}).(str_region).BLP = BLP_final;
        end
    end
end
disp('BLP analysis terminated ')


end

