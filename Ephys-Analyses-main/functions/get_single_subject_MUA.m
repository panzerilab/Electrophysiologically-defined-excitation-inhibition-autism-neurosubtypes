function MUA_data = get_single_subject_MUA(directory_name,subj_list,my_experiment_variables,subj,os_flag)
% This function gets single subject MUA data
% Inputs:
% directory name: the directory of the folder where the preprocessed files are 
% subj_list: cell array containing the name of the subjects
% my_experiment_variable: the structure defined in the ephys_analyses
% routine containing parameters of the experiment
% subj: the number of subject to read 
% os_flag: '/' for linux, '\' for windows
% Outputs: MUA_data: a structure containing the MUA data with the following
% fields:
% - fs: sampling frequency in Hz
% - type: type of the subject according to the subject file name (exp,ctrl,...)
% - MUA_data.(phase).(region): is a matrix of n_timepoints x n_channels 
% that contains the MUA data of the region (i.e. PFC) of the particular
% phase (i.e. baseline) according to the information in the
% my_experiment_variables 
% ---------------------------------------------------------------------
% this function reads the MUA of the subject #subj from the subj_list,
% chuncks the data into different phases, extracts the timing requested for
% each phase, and divides the data into all the channles based on the info
% in the my_experiment_variables that are coming from the excel file
% NOTE: BUG TO FIX: currently it is written in a way that has the
% assumption that the required part of the data is a multiplication of the
% length of the file. for example if the files are 5 minutes each and the
% user wants from minute 3 to 8, it does not work!!!
%% initialize 
MUA_data = struct;
subj_to_read = subj_list{subj};
path = strjoin({directory_name,subj_to_read,'MUA',''},os_flag);
file_list = dir(path);
n_files = length(file_list); % check how many files you have
phases = my_experiment_variables.phase;
n_phases = length(phases);
info = my_experiment_variables.data_info{subj};

for phase=1:n_phases
    all_data.(phases{phase}) = [];
end

%% detecting the phase for each file

for file=1:n_files % loop pver all files
    if ~file_list(file).isdir % only for files that are not directory
        temp_idx = file;
        file_name = file_list(file).name;
        for phase=1:n_phases % loop over phases
            if contains(file_name,phases{phase}) % check which phase does the file belong
               idx_file.(phases{phase})(file) = 1;
            end
        end
    end
end

%% determine the file length
file_name = file_list(temp_idx).name;
str_to_read = strcat(path,file_name);
temp = load(str_to_read);
fields = fieldnames(temp);
fs = temp.(fields{1}).fs;
len_data = size(temp.(fields{1}).data,1);
len_file = floor(len_data./fs./60);

for phase=1:n_phases % loop over phase and determine which files need to be read
    n_files_to_read.(phases{phase}) = info.(phases{phase}).duration/len_file;
end



%% reading files
for phase=1:n_phases
    tmp_idx = idx_file.(phases{phase});
    idx_exist = find(tmp_idx);
    tmp_n_files = n_files_to_read.(phases{phase});
    start_from = info.(phases{phase}).start;
    if start_from == -1 
        start_from = idx_exist(end) - tmp_n_files+1;
    else
        start_from = idx_exist(round(start_from/len_file)+1);
    end
    for file=1:tmp_n_files
        idx_to_read = start_from + file - 1;
        file_name = file_list(idx_to_read).name;
        str_to_read = strcat(path,file_name);
        temp = load(str_to_read);
        fields = fieldnames(temp);
        temp = temp.(fields{1}).data;
        all_data.(phases{phase}) = cat(1,all_data.(phases{phase}),temp); 
    end
end

%% extracting the chunck of required data,region assignment and bad channel removal
MUA_data.fs = fs;
regions = my_experiment_variables.regions;
n_regions = length(regions);
for phase=1:n_phases
    t = (info.(phases{phase}).duration*60)*fs;
    if size(all_data.(phases{phase}),1) >= t
        if info.(phases{phase}).start == 0
            t_start = size(all_data.(phases{phase}),1) - info.(phases{phase}).duration*60*fs;
            all_data.(phases{phase}) = all_data.(phases{phase})(t_start+1:end,:);
            MUA_data.time.(phases{phase}).ti = t_start/60;
            MUA_data.time.(phases{phase}).tf = (t_start + t)/60;
        else
            all_data.(phases{phase}) = all_data.(phases{phase})(1:t,:);
            MUA_data.time.(phases{phase}).ti = info.(phases{phase}).start;
            MUA_data.time.(phases{phase}).tf = info.(phases{phase}).start + info.(phases{phase}).duration;
        end
    end
    temp = all_data.(phases{phase});
    for reg=1:n_regions
        ch_to_select = my_experiment_variables.ch_list{subj,reg};
        temp_ch = temp(:,ch_to_select);
        ch_to_remove = my_experiment_variables.broken_ch_list{subj,reg};
        temp_ch(:,ch_to_remove) = [];
        MUA_data.(phases{phase}).(regions{reg}) = temp_ch;
    end
    
end

%% determining which group the subject belongs
groups = my_experiment_variables.groups;
group_names = my_experiment_variables.groups_name;
for grp=1:length(groups)
    if contains(file_name,groups{grp})
        MUA_data.type = group_names{grp};
    end  
end


end

