function [subj_list,data_info] = get_data_info(n_subjects,n_phases,phases,T)
% This funtion reads the excel containing the experiment info
% and extracts subject list and data lengths for further analyses
% Inputs: 
% n_subjects: number of subjects used for the analyses
% n_phases: number of phases you have 
% phases: cell array containing the name of the phases

subj_list = cell(n_subjects,1);
data_info = cell(n_subjects,1);

if n_phases ~= length(phases)
    error('number of phases does not match the phases list!!! fix it');
end

%T = readtable('Experiment_variables.xlsx');

for subj=1:n_subjects
    subj_list{subj} = T.subject(subj);
    tmp = struct;
    for ph=1:n_phases
        tmp.(phases{ph}).start = T.(strcat(phases{ph},'_start'))(subj);
        tmp.(phases{ph}).duration = T.(strcat(phases{ph},'_duration'))(subj);
    end
    data_info{subj} = tmp;
end

end

