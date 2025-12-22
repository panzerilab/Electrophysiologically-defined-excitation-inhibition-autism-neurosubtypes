function PLV_output = perform_PLV(LFP,plv_param,my_experiment_variables)
%   function to compute the phaselocking value (PLV) between channels of different regions
%   in the same phase
% Inputs:
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
% PLV_output : PLV_output.(phase_name).(str_region) is a structure that
% contains the PLV results for the phase = phase_name and region =
% region1_region2 (based on the my_experiment_variables, for example: PFC_rs) and contains the
% following fields:
%   - PLV: the calculated phase locking value between each channel in region 1 and each channel in region 2. the size
%         will be [1,n_channel_region1*n_channels_region2]
%   - phase_diff_mean: the angular mean of the difference of angles in time between each channel in region 1 and each channel in region 2. the size
%         will be [1,n_channels_region1*n_channels_region2]
%   - concentration: the concentration parameter of the angle differences between each channel in region 1 and each channel in region 2. the size
%         will be [1,n_channel_region1*n_channels_region2]
% --------------------------------------------------------------------
% at each phase and each region, the routine bandpass filters the signal using a
% non-causal user specified (filter order and bandwidth) butterwurth
% filter, then calculates the angle of the signal as the angle of 
% the Hilbert transformed signal. Then for each channel in region 1 and
% each channel in region2 it calculates the distribution of the angle
% differences by subtracting the angle of the region 2 from the region 1 at
% each timepoint, PLV is defined as the mean resultant vector (check circ_r) of this
% distribution, concentration is defined as the kappa (check circ_kappa)
% of the distribution, and the phase_diff_mean is the angular mean of the
% resulted angle differences distribution. 
% REQUIRES circular statistics toolbox
%% initialize
regions = my_experiment_variables.regions;
n_regions = length(regions);
phases = my_experiment_variables.phase;
n_phases = length(phases);
fs = LFP.fs;
freq = plv_param.frequency;
order = plv_param.filterorder;

if freq(1) == 0
    [b,a] = butter(order,freq(2)./(fs/2),'low');
else
    [b,a] = butter(order,freq./(fs/2),'bandpass');
end
PLV_output = struct;
%% calculate
for phase=1:n_phases % loop over phase
    for reg1=1:n_regions % loop over region 1 
        data_1 = LFP.(phases{phase}).(regions{reg1});
        
%         data_1_ch_avg = mean(data_1,2);
        data_1_LP = filtfilt(b,a,data_1); % filter signal in the required band
        analytical_sig_1 = hilbert(data_1_LP); % take the hilbert transform
        %angle_sig_1 = angle(analytical_sig_1); % calculate the angle
        
        for reg2=reg1+1:n_regions % loop over region 2 
            data_2 = LFP.(phases{phase}).(regions{reg2});
%             data_2_ch_avg = mean(data_1,2);
            data_2_LP = filtfilt(b,a,data_2); % filter signal in the required band
            analytical_sig_2 = hilbert(data_2_LP);  % take the hilbert transform
            %angle_sig_2 = angle(analytical_sig_2); % calculate the angle
            
            
            n_ch1 = size(analytical_sig_1,2);
            n_ch2 = size(analytical_sig_2,2);
            count = 1;
            phase_diff = zeros(n_ch1*n_ch2,length(analytical_sig_2)); % initialize the angle difference variable
            conc = zeros(n_ch1*n_ch2,1);  % initialize the variable to save the concentration parameter
            
            for ch1=1:n_ch1 % loop over channels in region 1 
                %data_ch1 = angle_sig_1(:,ch1); 
                for ch2=1:n_ch2 % loop over channels in region 2
                    %data_ch2 = angle(analytical_sig_2(:,ch2);
                    
                    phase_diff(count,:) = angle(analytical_sig_1(:,ch1)./analytical_sig_2(:,ch2));%data_ch1 - data_ch2; % take the angle differnce
                    conc(count) = circ_kappa(phase_diff(count,:)); % take the concentration parameter of angle differences
                    count = count + 1;
                end
            end
            region_1 = regions{reg1};
            region_2 = regions{reg2};
            str_region = strcat(region_1,'_',region_2);  % make the string as region1_region2 (i.e. PFC_rs)
            PLV_output.(phases{phase}).(str_region).PLV = circ_r(phase_diff'); % calculate the PLV value 
            PLV_output.(phases{phase}).(str_region).phase_diff_mean = circ_mean(phase_diff'); % calculate the angular mean of the angle differences
            PLV_output.(phases{phase}).(str_region).concentration = conc';
        end
    end  
end
disp('PLV analysis terminated ')



end

