function params = get_stat_params_LFP(list_of_statistics)
% function to extract the parameters for all the LFP stats requested in
% list_of_statistics
% Inputs:
% list_of_statistics: a cell array containing the name of the statistics
% requested
% Outputs:
% params: a structure with analyses names as its fields, each with their
% own subfields containing the parameters of the statistics requested, 
% for example: params.spectrum contains the %parameters for the spectrum
% statistics analyses
%% get parameters
n_statistics = length(list_of_statistics);
for stat=1:n_statistics
    if strcmp('spectrum',list_of_statistics{stat})
        params.spectrum.bands = {[0.1,1],[1,4],[4,8],[8,12],[12,30],[30,70]}; %requested frequency bands for statistics
        params.spectrum.bands_name = {'slow','delta','theta','alpha','beta','gamma'}; % corresponding bands name
        params.spectrum.fof = 1; % zero if you don't want the fof analyses
        params.spectrum.flim = 150; %Range of frequencies for cluster based correction
    elseif strcmp('coherency',list_of_statistics{stat})
        params.coherency.bands = {[0.1,1],[1,4],[4,8],[8,12],[12,30],[30,70]}; %requested frequency bands for statistics
        params.coherency.bands_name = {'slow','delta','theta','alpha','beta','gamma'}; % corresponding bands name; 
        params.coherency.flim = 40;   %Set frequency limits for cluster based analyses
    elseif strcmp('GC',list_of_statistics{stat})
        %params.GC.bands = {[0.1,1],[1,4],[4,8],[8,12],[12,30],[30,70]}; %requested frequency bands for statistics
        %params.GC.bands_name = {'infraslow','delta','theta','alpha','beta','gamma'}; % corresponding bands name;
        params.GC.bands = {[0.1,1],[1,4],[4,8],[8,12]};
        params.GC.bands_name = {'slow','delta','theta','alpha'};
    elseif strcmp('PLV',list_of_statistics{stat})
        params.PLV.sig_measure = {'PLV'}; % how you want to measure significancy of PLV, the other option is 'concentration';
        params.PLV.th = 0.1; % the threshold for assessing the significance {value}>th ==> significant
    elseif strcmp('slope',list_of_statistics{stat})
        params.slope.frequency = [20,40]; % the frequency you that you have fitted the line
    elseif strcmp('BLP',list_of_statistics{stat})
        params.BLP = [];
    else
        error('you have requested analyses which are not yet supported!!!');
    end
end


end

