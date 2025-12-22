function spikes = perform_spike_detection(MUA_data,params,my_experiment_variables)
% function to detect the spikes
% Inputs:
% MUA_data: a structure containing the LFP data with the following
% fields:
% - fs: sampling frequency in Hz
% - type: type of the subject according to the subject file name (exp,ctrl,...)
% - MUA_data.(phase_name).(region_name): is a matrix of n_timepoints x n_channels 
%   that contains the MUA data of the region (i.e. PFC) of the particular
%   phase (i.e. baseline) according to the information in the
%   my_experiment_variables

% params: a structure wich field names are the names of the required
%   analyses and inside each field the parameters of the analyses
% my_experiment_variables: the structure defined in the ephys_analyses
% routine containing parameters of the experiment
% Outputs:
% spikes : spikes.(phase_name).(str_region) is a structure that
% contains the spiking results for the phase = phase_name and region =
% region_name (based on the my_experiment_variables, for example: PFC) and contains the
% following fields:
%   - fs: the sampling rate for MUA
%   - times: a cell array of length n_channles, withing each the spike
%   times of that channel are stored in seconds
%   - idx: a cell array of length n_channles, withing each the spike idices
%   are stored
% --------------------------------------------------------------------
% this routine first loops over each region, concatenates the MUA from
% different phases, detects the threshold for spike detection based on the
% method specified in the params (check get_params_MUA.m) and
% (m_spikeDetector) for each region. Then loops over each phase and region
% and detects the spikes and outputs the results
% Note: Requires m_spikeDetector
%% initialize
regions = my_experiment_variables.regions;
phases = my_experiment_variables.phase;
n_regions = length(regions);
n_phases = length(phases);

method = params.spike_detect_method;
fs = MUA_data.fs;

spikes = struct;
th = struct;

%% set the same threshold for all phases 
for reg=1:n_regions % loop over region
    to_spike = [];
    for phase=1:n_phases
        to_spike = cat(1,to_spike,MUA_data.(phases{phase}).(regions{reg})); % concatanate the signal of different phases
    end
    [~, ~, ~, th.(regions{reg})] = m_spikeDetector(to_spike,fs,method); % detect the threshold for each region
end
clear to_spike;
%% detect spikes
for phase=1:n_phases % loop over phase
    for reg=1:n_regions % loop over region
        to_spike = MUA_data.(phases{phase}).(regions{reg});
        threshold = th.(regions{reg});
        [spk_times,spk_idx] = m_spikeDetector(to_spike,fs, '', threshold); % detect spikes
        spikes.(phases{phase}).(regions{reg}).times = spk_times;
        spikes.(phases{phase}).(regions{reg}).idx = spk_idx;
        spikes.fs = fs;
    end
end

end

