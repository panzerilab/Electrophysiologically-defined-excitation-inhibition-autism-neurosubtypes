function MUA_results = analyze_MUA(MUA_data,LFP_data,params,my_experiment_variables,subj)
%Function to perform multiple MUA analyses requested by the user
% Inputs:
% MUA_data: a structure containing the LFP data with the following
% fields:
% - fs: sampling frequency in Hz
% - type: type of the subject according to the subject file name (exp,ctrl,...)
% - MUA_data.(phase_name).(region_name): is a matrix of n_timepoints x n_channels 
%   that contains the MUA data of the region (i.e. PFC) of the particular
%   phase (i.e. baseline) according to the information in the my_experiment
%   variable
% LFP: a structure containing the LFP data with the following
% fields:
% - fs: sampling frequency in Hz
% - type: type of the subject according to the subject file name (exp,ctrl,...)
% - LFP_data.(phase_name).(region_name): is a matrix of n_timepoints x n_channels 
%   that contains the LFP data of the region (i.e. PFC) of the particular
%   phase (i.e. baseline) according to the information in the
%   my_experiment_variables
% params: a structure wich field names are the names of the required
%   analyses and inside each field the parameters of the analyses
% my_experiment_variables: the structure defined in the ephys_analyses
% routine containing parameters of the experiment
% Outputs:
% MUA_results :a structure with containing the results with each field representing 
%   the results of an analyse, for example: MUA_results.rate
%   contains the results of the rate
% Note: If PLV analyses is not required, the LFP can be passed as empty
% ([])
% Note: this function detects spikes from MUA before starting any analyses
%% initialize
MUA_results = struct;
analyze = params.analyses;
list_of_analyses = fieldnames(analyze); % get the name of the analyses
n_analyses = length(list_of_analyses); 

%% analyze
spikes = perform_spike_detection(MUA_data,params,my_experiment_variables); % detect spikes from MUA
for analysis=1:n_analyses % loop over each analyses
    analysis_to_do = list_of_analyses{analysis};
    if strcmp(analysis_to_do,'rate')
        MUA_results.(analysis_to_do) = perform_rate(spikes,MUA_data,analyze.(analysis_to_do),my_experiment_variables);
    elseif strcmp(analysis_to_do,'ISI')
        MUA_results.(analysis_to_do) = perform_ISI(spikes,my_experiment_variables);
    elseif strcmp(analysis_to_do,'PLV')
        MUA_results.(analysis_to_do) = perform_PLV_MUA(spikes,LFP_data,analyze.(analysis_to_do),my_experiment_variables,subj);
    else
        error('Requested analyses is not supported yet!!!')
    end
end

end

