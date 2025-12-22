function LFP_results = analyze_LFP(LFP,analyze,my_experiment_variables)
% Function to perform multiple LFP analyses
% Inputs:
% LFP: a structure containing the LFP data with the following
% fields:
% - fs: sampling frequency in Hz
% - type: type of the subject according to the subject file name (exp,ctrl,...)
% - LFP_data.(phase_name).(region_name): is a matrix of n_timepoints x n_channels 
%   that contains the LFP data of the region (i.e. PFC) of the particular
%   phase (i.e. baseline) according to the information in the my_experiment
%   variable
% analyze: a structure wich field names are the names of the required
%   analyses and inside each field the parameters of the analyses
% my_experiment_variable: the structure defined in the ephys_analyses
% routine containing parameters of the experiment
% Outputs:
% LFP_results :a structure with the following field: 
% - results: a structure containing the results with each field representing 
%   the results of an analyse, for example: LFP_results.results.coherency
%   contains the results of the coherency analyses


%% initialize
LFP_results = struct;

list_of_analyses = fieldnames(analyze);
n_analyses = length(list_of_analyses);
%% loop over analyses, check the ones that are requested, do them
for analysis=1:n_analyses
    analysis_to_do = list_of_analyses{analysis};
    if strcmp(analysis_to_do,'spectrum')
        LFP_results.(analysis_to_do) = perform_spectrum(LFP,analyze.(analysis_to_do),my_experiment_variables);
    elseif strcmp(analysis_to_do,'coherency')
        LFP_results.(analysis_to_do) = perform_coherency(LFP,analyze.(analysis_to_do),my_experiment_variables);
    elseif strcmp(analysis_to_do,'PLV')
        LFP_results.(analysis_to_do) = perform_PLV(LFP,analyze.(analysis_to_do),my_experiment_variables);
    elseif strcmp(analysis_to_do,'GC')
        LFP_results.(analysis_to_do) = perform_GC(LFP,analyze.(analysis_to_do),my_experiment_variables);
    elseif strcmp(analysis_to_do,'slope')
        LFP_results.(analysis_to_do) = perform_slope(LFP,analyze.(analysis_to_do),my_experiment_variables);
    elseif strcmp(analysis_to_do,'BLP')
        LFP_results.(analysis_to_do) = perform_BLP(LFP,analyze.(analysis_to_do),my_experiment_variables);
    else
        error('Requested analyses is not supported!!!')
    end
end


end


