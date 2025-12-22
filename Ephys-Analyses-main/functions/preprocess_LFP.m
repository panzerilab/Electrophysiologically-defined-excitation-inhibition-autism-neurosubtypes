function  LFP_out = preprocess_LFP(LFP_data,prprcss,my_experiment_variables)
% This function does the requested pre-processing on the LFP data
% Inputs:
% LFP_data: a structure containing the LFP data with the following
% fields:
% - fs: sampling frequency in Hz
% - type: type of the subject according to the subject file name (exp,ctrl,...)
% - LFP_data.(phase).(region): is a matrix of n_timepoints x n_channels 
%   that contains the LFP data of the region (i.e. PFC) of the particular
%   phase (i.e. baseline) according to the information in the my_experiment
%   variable
% - prprcss: a structure containing the parameters for the preprocesses,
%   look the function get_params_LFP for individual parameters for each preprocess 
% Outputs: 
%  - LFP_out: same as LFP_data after the preprocess has been applied  
% -------------------------------------------------------------------
% this function can perform 3 propecesses: 
% 1- adaptive line noise removal based on prep pipeline reported in:
% (https://www.frontiersin.org/articles/10.3389/fninf.2015.00016/full)
% 2- detrending the signal by dividing the signal into trials and removing
% the fitted line in the trial to remove trends in the signal that can be
% caused by movement based on prep pipeline reported in:
% (https://www.frontiersin.org/articles/10.3389/fninf.2015.00016/full)
% the parameters are:
%    preprocess.detrend.len_trial = 0.02;
% 3- Manual line noise removal by applying no phase shift notch filter using
% user specified BW and filter order on the user defined filter orders
% the parameters are:
%     prprcss.notch.order: a single number for the notch filter order
%     prprcss.preprocess.notch.f0: an array containing line frequencies to remove
%     prprcss.preprocess.notch.BW: a single number containing the bandwidth of the filter
% requires the EEGLAB and EEG-Clean-Tools-master
%% initialize parameters 
LFP_out = LFP_data;
list_of_preprocess = fieldnames(prprcss);
n_preprocess = length(list_of_preprocess);
phases = my_experiment_variables.phase;
regions = my_experiment_variables.regions;
n_phases = length(phases);
n_regions = length(regions);
fs = LFP_data.fs;
%% do preprocessing
for pr=1:n_preprocess % loop over preprocesses
    process_to_do = list_of_preprocess{pr};
    for phase=1:n_phases % loop over phase
        for reg=1:n_regions % loop over region
            Data_to_preprocess = LFP_out.(phases{phase}).(regions{reg});
            if strcmp('line',process_to_do)
                Data_to_line.data = Data_to_preprocess';
                Data_to_line.srate = fs;
                lineNoiseIn = prprcss.line;
                lineNoiseIn.Fs = fs;
                lineNoiseIn.lineNoiseChannels = 1:size(Data_to_preprocess,2); % apply line to all channels
                Data_out = cleanLineNoise(Data_to_line,lineNoiseIn);
                LFP_out.(phases{phase}).(regions{reg}) = Data_out.data';
            elseif strcmp('detrend',process_to_do) % detrend signal             
                det.detrendType = 'linear'; % linear detrending
                det.detrendChannels = 1:size(Data_to_preprocess,2); % all the channels
                det.detrendStepSize = prprcss.detrend.len_trial;
                Data_to_detrend.data = Data_to_preprocess';
                Data_to_detrend.srate = fs;
                Data_out = removeTrend(Data_to_detrend,det);
                LFP_out.(phases{phase}).(regions{reg}) = Data_out.data';
            elseif strcmp('notch',process_to_do)
                f0s = prprcss.notch.f0;
                Data_out = Data_to_preprocess;
                for freq=1:length(f0s)
                    f0 = f0s(freq);
                    BW = prprcss.notch.BW;
                    order = prprcss.notch.order;
                    notchspec = fdesign.notch('N,F0,BW',order,f0,BW,fs);
                    notchfilt = design(notchspec,'SystemObject',true);
                    Data_out = filtfilt(notchfilt.SOSMatrix,notchfilt.ScaleValues,Data_out);
                end
                LFP_out.(phases{phase}).(regions{reg}) = Data_out;
            else
                error('requested preprocess is not supported!!!');
            end
        end
    end
end

