function spectrum_output = perform_spectrum(LFP,spec,my_experiment_variables)
% Function to estimate power spectrum  
% Inputs:
% LFP: a structure containing the LFP data with the following
% fields:
% - fs: sampling frequency in Hz
% - type: type of the subject according to the subject file name (exp,ctrl,...)
% - LFP_data.(phase_name).(region_name): is a matrix of n_timepoints x n_channels 
%   that contains the LFP data of the region (i.e. PFC) of the particular
%   phase (i.e. baseline) 
% spec: a structure containing the parameters for the analyses, has the
% following fields:
%    - resolution: a positive number to determine the frequency resolution
%    - limits:  an array of lenght two ([f1,f2], f2>f1) for frequency limits to calculate
%    the spectrum
%    - pooling: 0 for no pooling, nan for automatic detection,other numbers for the minutes you want to pool
%    - ch_lvl:  0 for to average over channels, 1 for channel level analyses
% my_experiment_variable: the structure defined in the ephys_analyses
% routine containing parameters of the experiment
% Outputs:
% spectrum_output : spectrum_output.(phase_name).(region_name) is a structure that
% contains the spectrum results for the phase = phase_name and region =
% region_name (based on the my_experiment_variables) and contains the
% following fields:
%   - spectrum: the calculated raw power spectrum using pspectrum function. the size
%   will be [f1:freq_resoloution:f2] x (n_trials), n_trials is
%   determined based on pooling and wheter channel level analyses was 
%   requested or not. for example for a data of
%   10 minutes length and 32 channels, if the pooling is 1 and the ch_lvl
%   is 1, it will be 10*32, take note that first there will be the pooled versions 
%   of the same channel and then the next channel. if the ch_lvl = 0, it
%   will first calculate the power spectrum for each channel, then averages
%   over channels.
%   - spectrum_freq: an array of [f1:freq_resoloution:f2] x 1 which
%   contains the frequencies at which the spectrum was calculated

%% initialize
regions = my_experiment_variables.regions;
n_regions = length(regions);
phases = my_experiment_variables.phase;
n_phases = length(phases);
fs = LFP.fs;
freq_limits = spec.limits;
freq_res = spec.resolution;
%% calcualte
for reg=1:n_regions % loop over regions
    for phase=1:n_phases % loop over phases
        data = LFP.(phases{phase}).(regions{reg});
        data_avg = data;  
        n_ch = size(data_avg,2);
        if spec.pooling == 0 % check for pooling
            data_pool = data_avg;
        elseif isnan(spec.pooling)
            len_pool = automatic_pooling(data_avg,fs); 
            data_pool = reshape(data_avg,len_pool,[]);
        else
            len_pool = spec.pooling*60*fs;
            if rem(size(data_avg,1),len_pool) == 0
               data_pool = reshape(data_avg,len_pool,[]); 
            else
               error('the requested pooling length and data length are not compatible'); 
            end
        end
        n_trials = size(data_pool,2);
        fspec = [];
        Spec = [];
        for tr=1:n_trials
            LFP_to_Spectrogram = data_pool(:,tr);
            [Spec(:,tr),fspec(:,tr)] = pspectrum(LFP_to_Spectrogram, fs,'FrequencyLimits',freq_limits,'FrequencyResolution',freq_res);
            
        end
        if spec.ch_lvl == 0 % check to see for average over channel or not
           Spec = mean(reshape(Spec,size(Spec,1),[],n_ch),3);
        end
        spectrum_output.(phases{phase}).(regions{reg}).spectrum = Spec;
        spectrum_output.(phases{phase}).(regions{reg}).spectrum_freq = fspec(:,end);
    end
end
disp('Spectrum analysis terminated ')
end

