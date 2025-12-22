function ISI_output = perform_ISI(spikes,my_experiment_variables)
% function to compute the Inter spike intervals
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
% my_experiment_variable: the structure defined in the ephys_analyses
% routine containing parameters of the experiment
% Outputs:
% ISI_output : ISI_output.(phase_name).(str_region) is a structure that
% contains the ISI results for the phase = phase_name and region =
% region_name (based on the my_experiment_variables, for example: PFC) and contains the
% following fields:
%   - ISI: is a cell array of size 1 x n_channels which inside each cell is 
%     the inter spike interval in seconds for all the spikes of that
%     channel
%   - mean: is a cell array of length n_channles which inside each cell is
%       the average interspike interval for that channel 
% --------------------------------------------------------------------
% this routine, loops over each phase and each region, and then for each
% channel, calculates the difference of time in seconds between consecutive spikes. it
% stores both the distribution of the interspike intervals and their
% average for each channel seperately

%% initialize
regions = my_experiment_variables.regions;
n_regions = length(regions);
phases = my_experiment_variables.phase;
n_phases = length(phases);
ISI_output = struct;
%% calculate
for phase=1:n_phases % loop over phase
    for reg=1:n_regions % loop over region
        tmp_spikes = spikes.(phases{phase}).(regions{reg}).times;
        n_ch = length(tmp_spikes);
        for ch=1:n_ch % loop over channles
            tmp_ISI = diff(tmp_spikes{ch}); % calculate the differences between spike times 
            ISI_output.(phases{phase}).(regions{reg}).ISI{ch} = tmp_ISI; 
            ISI_output.(phases{phase}).(regions{reg}).mean(ch) = mean(tmp_ISI); % take the average ISI
        end
    end
end


end

