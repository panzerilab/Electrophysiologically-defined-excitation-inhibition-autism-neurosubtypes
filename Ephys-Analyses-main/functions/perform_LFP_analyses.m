function  perform_LFP_analyses(my_experiment_variables,options,directory_name,os_flag,pooling_flag)
%  This function performs the requested analyses on LFP data and saves the
%  results
% Inputs:
% my_experiment_variable: the structure defined in the ephys_analyses
% routine containing parameters of the experiment
% options: the structure defined in the ephys_analyses routine containing
%   requested analyses
% directory name: the directory of the preprocessed folder 
% os_flag: '/' for linux, '\' for windows
% Outputs:
% saves LFP_results seperately for each animal which is a structure with
% the following field: 
% - params: the parameters of the analyses (coming from get_params_LFP)
% - type: the subject type (exp or ctrl or ... ) based on filenames
% - results: a structure containing the results with each field representing 
%   the results of an analyse, for example: LFP_results.results.coherency
%   contains the results of the coherency analyses
%%                            
params = get_params_LFP(options); % get the parameters for all analyses
n_subjects = length(my_experiment_variables.subj_list);
LFP_results = struct;
LFP_results.params = params;


subj_name = dir(directory_name);
subj_list = {subj_name(3:end).name}; % get subjects

for subj=1:n_subjects % loop over each subject
    LFP_data = get_single_subject_LFP(directory_name,subj_list,my_experiment_variables,subj,os_flag); % read data
    LFP_results.type = LFP_data.type;
    LFP_results.time = LFP_data.time;  %Set the time proeprties of the analysis for the two phases
    if isfield(params,'preprocess') % if preprocess is required
        prprcss = params.preprocess;
        LFP_out = preprocess_LFP(LFP_data,prprcss,my_experiment_variables); % do the required preprocesses
    else
        LFP_out = LFP_data;
    end
    analyze = params.analyses;
    LFP_results.results = analyze_LFP(LFP_out,analyze,my_experiment_variables); % do the anayses
    save(strcat(pwd,os_flag,'results',os_flag,subj_list{subj},'_LFP_results.mat'),'LFP_results','-v7.3'); % save the results
    disp(['LFP analysis of subject',subj_list{subj},'terminated'])
end

end

