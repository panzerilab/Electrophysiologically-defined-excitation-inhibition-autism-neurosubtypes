function PLV_output = perform_PLV_MUA(spikes,LFP_data,plv_param,my_experiment_variables,subj)
% function to compute spike-LFP phase locking value (PLV)
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
% LFP: a structure containing the LFP data with the following
% fields:
% - fs: sampling frequency in Hz
% - type: type of the subject according to the subject file name (exp,ctrl,...)
% - LFP_data.(phase_name).(region_name): is a matrix of n_timepoints x n_channels 
%   that contains the LFP data of the region (i.e. PFC) of the particular
%   phase (i.e. baseline) 
% plv_param: a structure containing the parameters for the analyses, has the
% following fields:
%   - filterorder: integer determining the butterwurth filter order to use
%   - frequency = an array of lenght two of positive numbers ([f1,f2], f2>f1>=0) for the bandpass filter (if u want infraslow, use 0 instead 0.1 as the low frequency)
% my_experiment_variable: the structure defined in the ephys_analyses
% routine containing parameters of the experiment
% Outputs:
% PLV_output : PLV_output.(phase_name).(region_name) is a structure that
% contains the PLV results for the phase = phase_name and region =
% region_name (based on the my_experiment_variables, for example: PFC) and contains the
% following fields:
%   - PLV: the calculated phase locking value between spikes and the LFP of
%   the same channel. will be the length of n_channles
%   - phases: a cell array of the length of n_channels, inside each cell
%   is the distribution of the angles of the LFP of that channel when the same channel was spiking  
%   - mean: an array of 1 x n_channels, where the angular mean of the
%   distribution of LPF angles during spiking of each channel is stored 
% --------------------------------------------------------------------
% at each phase and each region, the routine bandpass filters the LPF using a
% non-causal user specified (filter order and bandwidth) butterwurth
% filter, then calculates the angle of the signal as the angle of 
% the Hilbert transformed signal. Then for each channel, it looks at the closest angle of 
% the LFP at time point of each spiketiming (using knn search) and creates the distribution of 
% LFP angles at spike times. the PLV is defined as the mean resultant
% vector of this distribution. the routine also outputs the distribtutions
% of the angels, and the circular average of the distribution
% REQUIRES circular statistics toolbox

%% initialize
regions = my_experiment_variables.regions;
n_regions = length(regions);
phases = my_experiment_variables.phase;
n_phases = length(phases);
fs_LFP = LFP_data.fs;
freq = plv_param.frequency;
order = plv_param.filterorder;

if freq(1) == 0
    [b,a] = butter(order,freq(2)./(fs_LFP/2),'low'); % define the filter
else
    [b,a] = butter(order,freq./(fs_LFP/2),'bandpass'); % define the filter
end
PLV_output = struct;

%% calculate
for phase=1:n_phases % loop over phase
    for reg=1:n_regions % loop over region
        tmp_spikes = spikes.(phases{phase}).(regions{reg}).times; % get the spikes
        data = LFP_data.(phases{phase}).(regions{reg}); % get the LFP
        data_LP = filtfilt(b,a,data); % filter the LFP
        analytical_sig = hilbert(data_LP); % hilbert transform the LFP 
        angle_sig = angle(analytical_sig); % get the angle
        len_LFP = size(angle_sig,1);
        time_LFP = linspace(0,len_LFP/fs_LFP,len_LFP); % generate the time for LFP
        n_ch = length(tmp_spikes);
        for ch=1:n_ch % loop over channels
            tmp_spikes_ch = tmp_spikes{ch};
            phase_tmp = angle_sig(:,ch); % get the channel angle
            phase_spike = phase_tmp(knnsearch(time_LFP',tmp_spikes_ch)); % find the angle of LFP at each spike time point
            PLV_output.(phases{phase}).(regions{reg}).phases{ch} = phase_spike;
            PLV_output.(phases{phase}).(regions{reg}).PLV(ch) = circ_r(phase_spike);
            PLV_output.(phases{phase}).(regions{reg}).mean(ch) = circ_mean(phase_spike);
        end
        %compute mean phase across channels
        %ch_median = circ_median(PLV_output.(phases{phase}).(regions{reg}).mean');
        %If the eman of one channel is more than 90 degrees different that
        %the mean across
        meanto2pi = wrapTo2Pi(PLV_output.(phases{phase}).(regions{reg}).mean);
        ch_to_flip = find(meanto2pi < pi/3 | meanto2pi > 5*pi/3);
       
        %ch_to_flip = find(abs(angdiff(PLV_output.(phases{phase}).(regions{reg}).mean,ch_median*ones(size(PLV_output.(phases{phase}).(regions{reg}).mean)))) > 2);
        %ch_to_flip = ch_to_flip(PLV_output.(phases{phase}).(regions{reg}).PLV(ch_to_flip) >0.1);
        PLV_output.(phases{phase}).(regions{reg}).ch_to_flip = ch_to_flip;
        disp(['Channels to flip in phase ',phases{phase},' in region ',regions{reg},'for subj ', ...
            my_experiment_variables.subj_list{subj}{1} ,':' ,num2str(ch_to_flip)]);
        PLV_output.(phases{phase}).(regions{reg}).channels = 1:n_ch;
    end
end


end

