function stats_all = perform_LFP_statistics(directory_name,os_flag)
% function to do the statistics on LFP results
% Inputs:
% directory name: the directory of the results folder 
% os_flag: '/' for linux, '\' for windows
% Outputs:
% stats_all: a structure with the name of the analyses as its fieldnames
% which contains the various group level statistical analyses on those
% results for example, stats_all.spectrum contains the statistics on the
% spectrum analyses
%-------------------------------------------------------------------------
% this routine loads the results of the LFP analyses (coming from ephys_analyses.m).
% puts them in a cell array of length of n_subjects. checks to see which
% analyses were done and does the group level statistics on all the
% analyses and outputs it in a structure named stats_all
%% initialize the list of subjects
subj_name = dir(directory_name);
subj_list = {subj_name(3:end).name};
idx_LFP = contains(subj_list,'LFP');
subj_list = subj_list(idx_LFP);

n_subjects = length(subj_list);
results_all = cell(1,n_subjects);
stats_all = struct;

%% read each subjects results
for subj=1:n_subjects
    path = strcat(directory_name,os_flag,subj_list{subj});
    results_all{subj} = load(path);
end

for subj=1:n_subjects
    results_all{subj} = results_all{subj}.LFP_results;
end

%% check the statistics needed and get the parameters
list_of_statistics = fieldnames(results_all{1}.params.analyses);  
n_statistics = length(list_of_statistics);
params = get_stat_params_LFP(list_of_statistics);

%% loop over each analyses and do the stats
for stat=1:n_statistics
    stat_to_do = list_of_statistics{stat};
    if strcmp(stat_to_do,'spectrum')
        stats_all.(stat_to_do) = perform_spectrum_stats(results_all,params.spectrum);
    elseif strcmp(stat_to_do,'coherency')
        stats_all.(stat_to_do) = perform_coherency_stats(results_all,params.coherency);
    elseif strcmp(stat_to_do,'PLV')
        stats_all.(stat_to_do) = perform_PLV_stats(results_all,params.PLV);
    elseif strcmp(stat_to_do,'GC')
        stats_all.(stat_to_do) = perform_GC_stats(results_all,params.GC);
    elseif strcmp(stat_to_do,'slope')
        stats_all.(stat_to_do) = perform_slope_stats(results_all,params.slope);
    elseif strcmp(stat_to_do,'BLP')
        stats_all.(stat_to_do) = perform_BLP_stats(results_all);
    end
end

