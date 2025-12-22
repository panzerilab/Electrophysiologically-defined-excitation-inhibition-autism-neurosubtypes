function [ch_list,broken_ch_list,toflip_ch_list] = get_channels(n_subjects,n_regions,regions,T,flip)
% This funtion reads the excel containing the experiment info
% and extracts channel list and broken channels for further analyses
% Inputs: 
% n_subjects: number of subjects used for the analyses
% n_regions: number of regions recorded
% regions: cell array containing the name of the name of the regions

ch_list = cell(n_subjects,n_regions);
broken_ch_list = cell(n_subjects,n_regions);
toflip_ch_list = cell(n_subjects,n_regions);

if n_regions ~= length(regions)
    error('number of regions does not match the regions list!!! fix it');
end
if nargin==3
    T = readtable('Experiment_variables.xlsx');
end

if length(T.subject) ~= n_subjects 
    error('number of subjects does not match the number subjects in the excel files!!! fix it')
end

last_ch_reg = zeros(n_subjects,n_regions);
for reg=1:n_regions
    reg_1 = T.Properties.VariableNames{reg+1};
    if ~strcmp(reg_1,regions{reg})
        error('there is a mismatch between region names in the files and in the script!!! fix it')
    else
        for subj=1:n_subjects
            
            str_br = strcat('br_',regions{reg});
            str_toflip = strcat('toflip_',regions{reg});
            temp_brch = T.(str_br)(subj);
            if ~strcmp(temp_brch,'~')
                broken_ch_list{subj,reg} = str2num(temp_brch{1});
            end
            if flip
                temp_toflip = T.(str_toflip)(subj);
                if ~strcmp(temp_toflip,'~')
                    toflip_ch_list{subj,reg} = str2num(temp_toflip{1});
                end
            end
            temp_ch_list = T.(regions{reg})(subj);
            if reg == 1
                ch_list{subj,reg} = last_ch_reg(subj,reg)+1:last_ch_reg(subj,reg)+temp_ch_list;
            else
                ch_list{subj,reg} = last_ch_reg(subj,reg-1)+1:last_ch_reg(subj,reg-1)+temp_ch_list;
            end
            last_ch_reg(subj,reg) = ch_list{subj,reg}(end);
        end
        
    end
end

end

