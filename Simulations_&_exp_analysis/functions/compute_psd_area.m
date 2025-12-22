function psd_area = compute_psd_area(signal,fs,f_range)
    % signal: Input signal
    % fs: Sampling frequency
    params.tapers=[3,5];
    [S,f] = mtspectrumc(signal,params);
    psd=S';

    fx=f*fs;
    idx=fx>f_range(1) & fx<f_range(2);

    fx=fx(idx);
    psd=psd(idx);
    % Take the logarithm of frequency and PSD

    % The slope of the line corresponds to the spectral slope
    psd_area = sum(psd);

end
