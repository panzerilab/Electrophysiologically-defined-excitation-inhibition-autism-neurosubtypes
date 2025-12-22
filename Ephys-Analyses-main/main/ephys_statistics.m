clear all
clc
close all

%% add directory
os_flag = filesep; % file seperator of the operating system 
mydir = pwd;
idcs   = strfind(mydir,os_flag);
newdir = mydir(1:idcs(end)-1);
addpath(genpath(newdir));

%% set options
options.LFP.do = 1; %1 if u want the statistics on LFP
options.MUA.do = 0;%1 if u want the statistics on MUA
directory_name = uigetdir(pwd, 'Select results folder'); % where the results folder is saved
%% Do statistics
if options.LFP.do
    stats_LFP = perform_LFP_statistics(directory_name,os_flag);
end

if options.MUA.do
    rt.ch_lvl = 1; % 1 to do ch_lvl, 0 for channel average
    rt.pooling = 0; %  for requested minutes, 0 for no pooling
    rt.eff = 'higher'; % 'lower' if u expect rate going low, 'higher' if u expect it to go high
    stats_MUA = perform_MUA_statistics(directory_name,rt,os_flag);
end