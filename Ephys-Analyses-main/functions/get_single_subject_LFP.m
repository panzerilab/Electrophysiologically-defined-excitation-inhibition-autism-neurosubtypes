function LFP_data = get_single_subject_LFP(directory_name,subj_list,my_experiment_variables,subj,os_flag)
% This function gets single subject LFP data
% Inputs:
% directory name: the directory of the preprocessed folder 
% subj_list: cell array containing the name of the subjects
% my_experiment_variable: the structure defined in the ephys_analyses
% routine containing parameters of the experiment
% subj: the number of subject to read 
% os_flag: '/' for linux, '\' for windows
% Outputs: LFP_data: a structure containing the LFP data with the following
% fields:
% - fs: sampling frequency in Hz
% - type: type of the subject according to the subject file name (exp,ctrl,...)
% - LFP_data.(phase).(region): is a matrix of n_timepoints x n_channels 
% that contains the LFP data of the region (i.e. PFC) of the particular
% phase (i.e. baseline) according to the information in the
% my_experiment_variables 
% ---------------------------------------------------------------------
% this function reads the LFP of the subject #subj from the subj_list,
% chuncks the data into different phases, extracts the timing requested for
% each phase, and divides the data into all the channles based on the info
% in the my_experiment_variables that are coming from the excel file

%% initialize
LFP_data = struct;
subj_to_read = subj_list{subj};
path = strjoin({directory_name,subj_to_read,'LFP',''},os_flag);
file_list = dir(path);
n_files = length(file_list);
phases = my_experiment_variables.phase;
n_phases = length(phases);
info = my_experiment_variables.data_info{subj};

for phase=1:n_phases
    all_data.(phases{phase}) = [];
end

%% reading data
for file=1:n_files
    if ~file_list(file).isdir
        file_name = file_list(file).name;
        str_to_read = strcat(path,file_name);
        temp = load(str_to_read);
        fields = fieldnames(temp);
        fs = temp.(fields{1}).fs;
        if strcmp(my_experiment_variables.data_to_read,'No50Hz')
            temp = temp.(fields{1}).dataNo50Hz;
        else
            temp = temp.(fields{1}).data;
        end
        for phase=1:n_phases
            if contains(file_name,phases{phase})
               all_data.(phases{phase}) = cat(1,all_data.(phases{phase}),temp); 
            end
        end
        
    end
    
end

%% extracting the chunck of required data,region assignment and bad channel removal
LFP_data.fs = fs;
regions = my_experiment_variables.regions;
n_regions = length(regions);
for phase=1:n_phases
    t = (0:1/fs:(length(all_data.(phases{phase}))-1)/fs);
    disp(['Total duration of phase ', phases{phase},':',num2str(t(end)/60)])
    info_phase = info.(phases{phase});
    %Check for special case below: If start is set at "-1" use 
    %the last n-minutes, not n-minutes from the starting time
    if info_phase.start == -1
        t_end = t(end);
        t_start = t(end) - info_phase.duration*60;
        LFP_data.time.(phases{phase}).ti = t_start/60;
        LFP_data.time.(phases{phase}).tf = t(end)/60;
    else
        t_start = (info_phase.start)*60;
        t_end = t_start + info_phase.duration*60;
        LFP_data.time.(phases{phase}).ti = info_phase.start;
        LFP_data.time.(phases{phase}).tf = info_phase.start + info_phase.duration;   
    end
    disp(['Analyze from minute ', num2str(t_start/60), ' to minute ',num2str(t_end/60)])
    idx_to_extract = t<t_end & t>=t_start;
    if rem(sum(idx_to_extract),10) == 9
        tmp_idx_to_fix = find(idx_to_extract,1,'last');
        if tmp_idx_to_fix + 1 <= length(idx_to_extract)
            idx_to_extract(tmp_idx_to_fix +1) = 1;
        else
            tmp_idx_to_fix = find(idx_to_extract,1,'first');
            idx_to_extract(tmp_idx_to_fix - 1) = 1;
        end
    elseif rem(sum(idx_to_extract),10) == 1
        tmp_idx_to_fix = find(idx_to_extract,1,'last');
        idx_to_extract(tmp_idx_to_fix) = 0;
    end
    temp = all_data.(phases{phase});
    temp = temp(idx_to_extract,:);
    for reg=1:n_regions % select channels
        ch_to_select = my_experiment_variables.ch_list{subj,reg};
        temp_ch = temp(:,ch_to_select);
        ch_to_remove = my_experiment_variables.broken_ch_list{subj,reg};
        temp_ch(:,ch_to_remove) = []; % remove bad channels
        if my_experiment_variables.flip ==1
            %Flip channel on the basis of the excel sheet
            ch_to_flip = my_experiment_variables.toflip_ch_list{subj,reg};
            temp_ch(:,ch_to_flip) = - temp_ch(:,ch_to_flip);
        end
        LFP_data.(phases{phase}).(regions{reg}) = temp_ch;
    end
end

%% determining which group the subject belongs
groups = my_experiment_variables.groups;
group_names = my_experiment_variables.groups_name;
for grp=1:length(groups)
    if contains(file_name,groups{grp})
        LFP_data.type = group_names{grp};
    end  
end


end

