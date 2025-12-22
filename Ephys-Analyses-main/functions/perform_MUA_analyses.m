function perform_MUA_analyses(my_experiment_variables,options,directory_name,os_flag)
%  This function performs the requested analyses on MUA data one subject at a time and saves the
%  results
% Inputs:
% my_experiment_variable: the structure defined in the ephys_analyses
% routine containing parameters of the experiment
% options: the structure defined in the ephys_analyses routine containing
%   requested analyses
% directory name: the directory of the preprocessed folder 
% os_flag: '/' for linux, '\' for windows
% Outputs:
% saves MUA_results seperately for each animal which is a structure with
% the following field: 
% - params: the parameters of the analyses (coming from get_params_MUA)
% - type: the subject type (exp or ctrl or ... ) based on filenames
% - results: a structure containing the results with each field representing 
%   the results of an analyse, for example: MUA_results.results.rate
%   contains the results of the rate analyses

%%
params = get_params_MUA(options); % get the parameters for all MUA analyses
n_subjects = length(my_experiment_variables.subj_list);
MUA_results = struct;
MUA_results.params = params;


subj_name = dir(directory_name);
subj_list = {subj_name(3:end).name}; % get subjects

for subj=1:n_subjects % loop over subjects
    MUA_data = get_single_subject_MUA(directory_name,subj_list,my_experiment_variables,subj,os_flag);
    if sum(contains(options.MUA.analyses,'PLV'))>0 % check if PLV analyses is requested, so it will need the LFP
       LFP_data = get_single_subject_LFP(directory_name,subj_list,my_experiment_variables,subj,os_flag); % read LFP
       if options.LFP.preprocess.do == 1 % chack if the LFP needs preprocessing
           params_LFP = get_params_LFP(options); % get parameters for LFP preprocessing
           prprcss = params_LFP.preprocess;
           LFP_data = preprocess_LFP(LFP_data,prprcss,my_experiment_variables); % preprocess the LFP
       end
    else
        LFP_data = [];
    end
    
    MUA_results.type = MUA_data.type; % set the type (exp,ctrl,..) of the data
    MUA_results.time = MUA_data.time; %Set the time proeprties of the analysis for the two phases
    MUA_results.results = analyze_MUA(MUA_data,LFP_data,params,my_experiment_variables,subj); % analyze MUA data
    save(strcat(pwd,os_flag,'results',os_flag,subj_list{subj},'_MUA_results.mat'),'MUA_results','-v7.3'); % save results
    disp(['MUA analysis of subject',subj_list{subj},'terminated'])
end


end

