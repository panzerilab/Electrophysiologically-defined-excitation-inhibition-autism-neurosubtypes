function stats_all = perform_MUA_statistics(directory_name,rt,os_flag)
% function to do the statistics on MUA results
% Inputs:
% directory name: the directory of the results folder 
% rt: a structure containing some details for rate analyses with the
% following fields:
%  - ch_lvl: 1 to do ch_lvl analyses, 0 for channel average 
%  - pooling an integer value to pool data in minutes, 0 for no pooling
%  - eff: 'lower' if u expect rate going low, 'higher' if u expect it to go high
% os_flag: '/' for linux, '\' for windows
% Outputs:
% stats_all: a structure with the name of the analyses as its fieldnames
% which contains the various group level statistical analyses on those
% results for example, stats_all.rate contains the statistics on the
% rate analyses
%-------------------------------------------------------------------------
% this routine loads the results of the LFP analyses (coming from ephys_analyses.m).
% puts them in a cell array of length of n_subjects. checks to see which
% analyses were done and does the group level statistics on all the
% analyses and outputs it in a structure named stats_all


%% get subjects info
subj_name = dir(directory_name);
subj_list = {subj_name(3:end).name};
idx_MUA = contains(subj_list,'MUA');
subj_list = subj_list(idx_MUA);

n_subjects = length(subj_list);
results_all = cell(1,n_subjects);

stats_all = struct;
%% read results
for subj=1:n_subjects
    path = strcat(directory_name,os_flag,subj_list{subj});
    results_all{subj} = load(path);
end

list_of_statistics = fieldnames(results_all{1}.MUA_results.params.analyses); % get the name of required statistics
n_statistics = length(list_of_statistics);

%% do stats
for stat=1:n_statistics % loop over statistics
    stat_to_do = list_of_statistics{stat};
    if strcmp(stat_to_do,'rate')
        stats_all.(stat_to_do) = perform_rate_stats(results_all,rt);
    elseif strcmp(stat_to_do,'ISI')
        stats_all.(stat_to_do) = perform_ISI_stats(results_all);
    elseif strcmp(stat_to_do,'PLV')
        stats_all.(stat_to_do) = perform_PLV_MUA_stats(results_all);
    end
end


