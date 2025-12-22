function slope_output = perform_slope(LFP,slope_param,my_experiment_variables)
%   function to compute the slope of spectrum
% Inputs:
% LFP: a structure containing the LFP data with the following
% fields:
% - fs: sampling frequency in Hz
% - type: type of the subject according to the subject file name (exp,ctrl,...)
% - LFP_data.(phase_name).(region_name): is a matrix of n_timepoints x n_channels 
%   that contains the LFP data of the region (i.e. PFC) of the particular
%   phase (i.e. baseline) 
% slope_param: a structure containing the parameters for the analyses, has the
% following fields:
%    - resolution: a positive number to determine the frequency resolution
%    - limits:  an array of lenght two ([f1,f2], f2>f1) for frequency limits to calculate
%    the spectrum
%    - pooling: 0 for no pooling, nan for automatic detection,other numbers for the minutes you want to pool
%    - ch_lvl:  0 for to average over channels, 1 for channel level analyses
%    - frequency: an array of lenght two ([f1,f2], f2>f1) for frequency
%       limits to fit the line
% my_experiment_variable: the structure defined in the ephys_analyses
% routine containing parameters of the experiment
% Outputs:
% slope_output : slope_output.(phase_name).(region_name) is a structure that
% contains the slope results for the phase = phase_name and region =
% region_name (based on the my_experiment_variables) and contains the
% following fields:
%   - spectrum: the calculated raw power spectrum using pspectrum function. the size
%   will be [f1:freq_resoloution:f2] x (n_trials), n_trials is
%   determined based on pooling and whether channel level analyses was 
%   requested or not. for example for a data of
%   10 minutes length and 32 channels, if the pooling is 1 and the ch_lvl
%   is 1, it will be 10*32, take note that first there will be the pooled versions 
%   of the same channel and then the next channel. if the ch_lvl = 0, it
%   will first calculate the power spectrum for each channel, then averages
%   over channels.
%   - spectrum_freq: an array of [f1:freq_resoloution:f2] x 1 which
%   contains the frequencies at which the spectrum was calculated
%   - slope: a cell array of size 1 x n_trials which inside each cell is an array of 
%     of length 2 containing the slope and bias of the fitted line
% --------------------------------------------------------------------
% at each phase and each region, the routine calculates the spectrum using
% perform_spectrum (read the documentation) and fits a line in the user
% specified frequency by having the log10 of power as y and log10 of the
% frequency as x
% REQUIRES perform_spectrum.m 

%% initialize
regions = my_experiment_variables.regions;
n_regions = length(regions);
phases = my_experiment_variables.phase;
n_phases = length(phases);
freq_to_fit = slope_param.frequency;
slope_output = struct;

%% calculate
spectrum_output = perform_spectrum(LFP,slope_param,my_experiment_variables); % calculate spectrum
for reg=1:n_regions % loop over regions 
    for phase=1:n_phases % loop over phases
        data_pool = spectrum_output.(phases{phase}).(regions{reg}).spectrum;
        f_spec = spectrum_output.(phases{phase}).(regions{reg}).spectrum_freq;
        n_trials = size(data_pool,2);
        upp_lim = freq_to_fit(2);
        low_lim = freq_to_fit(1);
        beta = cell(1,n_trials); % initialize the line coefficients
        for tr=1:n_trials % loop over trials
            S = data_pool(:,tr);
            f = f_spec;
            idx_line_fit = f < upp_lim & f > low_lim; % find the indices to fit the line
            x = log10(f(idx_line_fit)); % go to log10
            y = log10(S(idx_line_fit));% go to log10
            beta{tr} = polyfit(x,y,1); % fit the line
        end
        slope_output.(phases{phase}).(regions{reg}).spectrum = data_pool;
        slope_output.(phases{phase}).(regions{reg}).spectrum_freq = f_spec;
        slope_output.(phases{phase}).(regions{reg}).slope = beta;
    end
end

disp('Slope analysis terminated ')


end

