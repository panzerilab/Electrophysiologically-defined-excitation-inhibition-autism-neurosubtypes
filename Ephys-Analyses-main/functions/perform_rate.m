function rate_output = perform_rate(spikes,MUA_data,rt,my_experiment_variables)
% function to compute the rate properties
% Inputs:
% spikes : spikes.(phase_name).(str_region) is a structure that
% contains the spiking results for the phase = phase_name and region =
% region_name (based on the my_experiment_variables, for example: PFC) and contains the
% following fields:
%   - fs: the sampling rate for MUA
%   - times: a cell array of length n_channles, withing each the spike
%   times of that channel are stored in seconds
%   - idx: a cell array of length n_channles, withing each the spike idices
%   are stored
% MUA_data: a structure containing the LFP data with the following
% fields:
% - fs: sampling frequency in Hz
% - type: type of the subject according to the subject file name (exp,ctrl,...)
% - MUA_data.(phase_name).(region_name): is a matrix of n_timepoints x n_channels 
%   that contains the MUA data of the region (i.e. PFC) of the particular
%   phase (i.e. baseline) according to the information in the
%   my_experiment_variables

% rt: a structure containing the parameters for the analyses, has the
% following fields:
%       - time_bin: integer value in seconds to bin the data and calculate the rate
% my_experiment_variable: the structure defined in the ephys_analyses
% routine containing parameters of the experiment
% Outputs:
% rate_output : rate_output.(phase_name).(str_region) is a structure that
% contains the rate results for the phase = phase_name and region =
% region_name (based on the my_experiment_variables, for example: PFC) and contains the
% following fields:
%   - rate: the calculated spike rate (in n_spikes/seconds) will be the size of n_bins x
%   n_channels. n_bins depends on the length of data and the requested time
%   bin. for example if the data is 10 minutes and the requested timebin is
%   60 seconds, it will be 10
%   - time_bin: the time_bin used to bin the data and calculating the rate
% --------------------------------------------------------------------
% this routine, loops over each phase and each region, and then for each
% channel, segments the data in time into equel timebins (according to the
% requested timebin), finds the number of spikes in each segments, divides
% it by the timebin to calculate the rate at each timebin
   

%% initialize
regions = my_experiment_variables.regions;
n_regions = length(regions);
phases = my_experiment_variables.phase;
n_phases = length(phases);
fs = spikes.fs;
time_bin = rt.time_bin;
rate_output = struct;

%% calculate 
for phase=1:n_phases % loop over phase
    for reg=1:n_regions % loop pver region
        spk = spikes.(phases{phase}).(regions{reg}).idx; % get the indices of spikes
        spike_idx = zeros(size(MUA_data.(phases{phase}).(regions{reg})));
        [len_data,n_ch] = size(spike_idx);
        len_bin = time_bin*fs;
        if rem(len_data,len_bin) > 0 % check the data length and binning rate are compatible
            error('data length and binning length are not compatible')
        else
            len_rate = len_data/len_bin;
            rate = zeros(len_rate,n_ch); % create the rate matrix
            for ch=1:n_ch % loop over channels
                count = spk{ch};
                spike_idx(count,ch) = 1;
                tmp_count = reshape(spike_idx(:,ch),[],len_rate); % reshape based on requested timebin
                rate(:,ch) = sum(tmp_count)./time_bin; % calculate rate
            end
        end
        rate_output.(phases{phase}).(regions{reg}).rate = rate;
        rate_output.(phases{phase}).(regions{reg}).time_bin = time_bin;
    end
end



end

