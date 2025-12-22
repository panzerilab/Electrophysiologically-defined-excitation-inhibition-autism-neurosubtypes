function params = get_params_LFP(options,pooling)
% function to extract the parameters for all the LFP analyses requested in
% options
% Inputs:
% options: the structure defined in the ephys_analyses
% Outputs:
% params: a structure with fields preprocess and analyses, each with their
% own subfields containing the name of the preprocess or the analyses and
% its parameters, for example: params.analyses.spectrum contains the
% parameters for the spectrum analyses


%% pre-processing
params = struct;
if options.LFP.preprocess.do
    list_of_preprocess = options.LFP.preprocess.analyses;
    n_preprocess = length(list_of_preprocess);
    for prcss=1:n_preprocess
        if strcmp('line',list_of_preprocess{prcss})
             params.preprocess.line = struct('lineFrequencies', [50,100],...
                'taperBandWidth',2,'taperWindowSize',4,'taperWindowStep',2,'tau',100,...
                'maximumIterations',100,'fPassBand',[40,120],'fScanBandWidth',10,'p',0.05,'pad',0);
            
        elseif strcmp('detrend',list_of_preprocess{prcss})
            params.preprocess.detrend.len_trial = 0.02; % seconds for dividing into trials to apply detrending 
        elseif strcmp('notch',list_of_preprocess{prcss})
            params.preprocess.notch.order = 2; % notch filter order
            params.preprocess.notch.f0 = [50,100]; % line frequencies to remove
            params.preprocess.notch.BW = 0.4; % the bandwidth of the filter
            
        else 
            error('you have requested a preprocess which is not supported!!!');
        end
    end
end
%% analyses
list_of_analyses = options.LFP.analyses;
n_analyses = length(list_of_analyses);
for analys=1:n_analyses
    if strcmp('spectrum',list_of_analyses{analys})
        params.analyses.spectrum.resolution = 1; % required frequency resolution
        params.analyses.spectrum.limits = [0,250]; % frequency limits to calculate
        params.analyses.spectrum.pooling = nan; % 0 for no pooling, nan for automatic detection,other numbers for the minutes you want to pool
        params.analyses.spectrum.ch_lvl = 0; % 0 for channel average, 1 for channel level analyses
    elseif strcmp('coherency',list_of_analyses{analys})
        params.analyses.coherency.resolution = 0.1; % required frequency resolution
        params.analyses.coherency.limits = [0,100]; % frequency limits to calculate
        params.analyses.coherency.pooling = nan; % zeros for no pooling, nan for automatic detection,other numbers for the minutes you want to pool
        params.analyses.coherency.window = 2000; % integer value used for windowing the data to estimate the coherency 
        params.analyses.coherency.noverlap = 1000; % integer value used for amount of overlapping between windows (must be smaller than window)
        params.analyses.coherency.ch_lvl = 0; % 0 for channel average, 1 for channel level analyses
        params.analyses.coherency.confidence = 0; % 1 for extimating the coherence confidence level with bootstrap methods (default 0.95)
    elseif strcmp('GC',list_of_analyses{analys})
        params.analyses.GC.down_sampling_factor = 25;
        params.analyses.GC.pooling = 0; % zeros for no pooling,other numbers for the minutes you want to pool (automatic pooling is not supported for cGC)
        params.analyses.GC.len_trial = 1; % minutes for dividing into trials 
        params.analyses.GC.resolution = 0.1; % required frequency resolution
        params.analyses.GC.icregmode = 'LWR'; % type of regression for estimating order
        params.analyses.GC.regmode = 'OLS'; % type of regression for final GC computing
        params.analyses.GC.alpha = 0.05; % significance level for testing causality in time
        params.analyses.GC.mhtc = 'FDR';% multiple correction method for CG in time 
        params.analyses.GC.n_shuffles = 50; % number of shuffles to estimate GC in frequency
        params.analyses.GC.max_order = 15; % maximum order to consider for AR model
        params.analyses.GC.ch_lvl = 1;
    elseif strcmp('PLV',list_of_analyses{analys})
        params.analyses.PLV.filterorder = 3; % filter order to use
        params.analyses.PLV.frequency = [0,2]; % frequency to consider for the bandpass filtering (if u want slow, use 0 instead 0.1 as the low frequency)
    elseif strcmp('slope',list_of_analyses{analys})
        params.analyses.slope.frequency = [20,40]; % the frequency you want to fit the line
        params.analyses.slope.resolution = 2; % required frequency resolution 
        params.analyses.slope.limits = [0,100]; % frequency limits to calculate
        params.analyses.slope.pooling = nan; % 0 for no pooling, nan for automatic detection,other numbers for the minutes you want to pool
        params.analyses.slope.ch_lvl = 0; % 0 for channel average, 1 for channel level analyses
    elseif strcmp('BLP',list_of_analyses{analys})
        params.analyses.BLP.frequency = [30,70]; % the frequency you want the BLP to be calculated
        params.analyses.BLP.filterorder = 3; % filter order to use
        params.analyses.BLP.pooling = 10; % 0 for no pooling, nan for automatic detection,other numbers for the minutes you want to pool
        %BLP calculate coherence in a slow band with limits 0-1 Hz and
        %parameters as follows
        params.analyses.BLP.window = 2000;   
        params.analyses.BLP.noverlap = 1000;
        params.analyses.BLP.limits = [0.01,1];
        params.analyses.BLP.resolution = 0.004;
    else
        error('you have requested analyses which are not yet supported!!!');
    end
end


end

