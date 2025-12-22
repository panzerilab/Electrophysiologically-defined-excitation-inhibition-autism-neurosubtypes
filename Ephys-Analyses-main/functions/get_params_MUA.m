function params = get_params_MUA(options)
% function to extract the parameters for all the MUA analyses requested in
% options
% Inputs:
% options: the structure defined in the ephys_analyses
% Outputs:
% params: a structure with fields preprocess and analyses, each with their
% own subfields containing the name of the preprocess or the analyses and
% its parameters, for example: params.analyses.rate contains the
% parameters for the rate analyses

params = struct;
%% pre-processing

%% analyses
params.spike_detect_method = 'median';  % method for spike detection, acceptable options are "median", "5std", 8std
list_of_analyses = options.MUA.analyses;
n_analyses = length(list_of_analyses);
for analys=1:n_analyses
    if strcmp('rate',list_of_analyses{analys})
        params.analyses.rate.time_bin = 60; % time bin to calculate the rate in seconds
    elseif strcmp('ISI',list_of_analyses{analys})
        params.analyses.ISI = []; 
    elseif strcmp('PLV',list_of_analyses{analys})
        params.analyses.PLV.filterorder = 3; % filter order to use
        params.analyses.PLV.frequency = [0,2.4]; % frequency to consider for the bandpass filtering (if u want infraslow, use 0 instead 0.1 as the low frequency)
    else
        error('you have requested analyses which are not supported yet!!!');
    end
end


end