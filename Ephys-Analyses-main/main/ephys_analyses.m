clear all
clc
close all

%% adding paths
mydir = pwd;
idcs   = strfind(mydir,filesep);
newdir = mydir(1:idcs(end)-1);
addpath(genpath(newdir)); % adding all the necessary folders for all the analyses
set_up_toolboxes();
clear

mkdir results
os_flag = filesep;
%% pre defining vriables (variables related to the experiment)
my_experiment_variables = struct; % this structure contains variables related to the experiment
my_experiment_variables.groups = {'exp','ctrl'}; % here you put the names of the different groups you have according to the "FILENAMES"!!! 
my_experiment_variables.groups = {'exp','ctrl'}; % here you put the names of the different groups you have according to the "FILENAMES"!!! 
my_experiment_variables.groups_name = {'exp','ctrl'}; % here you put the names of groups based on what you want to be shown on figures
my_experiment_variables.regions = {'PFC','rs'}; % here you put the names of regions based on what you want to be shown on figures
%my_experiment_variables.phase = {'baseline','CNO'}; % here you put the names of the different phases you have according to the "FILENAMES"!!!
%%
%Load all metadata (then choose which one to analyse)
TCAMK = readtable('CAMKExperiment_variables_flip=0-2.xlsx');
TPV = readtable('PVExperiment_variables_flip=0-2.xlsx');
THSYN = readtable('Experimental_variables_3electrode.xlsx');
%TCAMK_s = readtable('Camk_round2_Experiment_variables.xlsx');
%TPV_new = readtable('PVExperiment_variables_newctrl.xlsx');
myT = TCAMK;
%phase1 = split(myT.Properties.VariableNames{6},'_');
%phase2 = split(myT.Properties.VariableNames{7},'_');
phase1 = split(myT.Properties.VariableNames{8},'_');
phase2 = split(myT.Properties.VariableNames{9},'_');
my_experiment_variables.phase = {phase1{1},phase2{1}};
%% assiging channels and broken channels etc
n_subjects =  size(myT,1); % number of subjects, put one if you are just checking a subject
n_regions = length(my_experiment_variables.regions); % calculates how many regions you have
my_experiment_variables.flip = 0; %D if flip channels based on the analyses of spikes locking value
% extracts the channel list and broken channel lists  and eventually the
% channels to flip
[my_experiment_variables.ch_list,my_experiment_variables.broken_ch_list,my_experiment_variables.toflip_ch_list] = get_channels(n_subjects,n_regions,my_experiment_variables.regions,myT,my_experiment_variables.flip); 
%% defining variables related to reading data
n_phases = length(my_experiment_variables.phase);
[my_experiment_variables.subj_list,my_experiment_variables.data_info] = get_data_info(n_subjects,n_phases,my_experiment_variables.phase,myT); % extracts the timing information
my_experiment_variables.data_to_read = 'No50Hz';%'raw'; % the other option is 'No50Hz' which gives you the no 50Hz data
 %Decide whethere to flip time series on the basis of PLV MUA analyses
%% defining options for what type of analyses to perform
options = struct; % this structure contains variables related to the analyses

options.LFP.do = 1; % 0 if u don't want the LFP to be analyzed
options.MUA.do = 0; % 0 if u don't want the MUA to be analyzed

options.LFP.analyses = {'spectrum','coherency','PLV','GC','slope','BLP'}; % list of the analyses you want to do with LFP
options.MUA.analyses = {'rate','ISI','PLV'}; % list of the analyses you want to do with MUA

options.LFP.preprocess.do = 0; % if u want to preprocess LFP, else put 0; (still needs to be refined)
options.LFP.preprocess.analyses = {'line','detrend','notch'}; % list of preprocessing u want to do on LFP


%% Do analyses

directory_name = uigetdir(pwd, 'Select experiment preprocessed folder');
if options.LFP.do
   perform_LFP_analyses(my_experiment_variables,options,directory_name,os_flag); 
end

if options.MUA.do
   perform_MUA_analyses(my_experiment_variables,options,directory_name,os_flag);
end
