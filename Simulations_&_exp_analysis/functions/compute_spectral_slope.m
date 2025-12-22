function spectral_slope = compute_spectral_slope(signals, fs, f_range)
    % signals: Matrix of input signals (each signal is a column)
    % fs: Sampling frequency
    % f_range: Frequency range for analysis

    num_signals = size(signals, 2);
    spectral_slope = zeros(1, num_signals);
    params.tapers=[3,5];

    for i = 1:num_signals
        signal = signals(:, i);

        % Compute spectral slope for each signal
        [S, f] = mtspectrumc(signal,params);
        psd = S';
        fx = f * fs;
        
        % Interpolate to avoid drop around 55 Hz (European powerline)
        drop_band = [44, 56];  % Hz range to interpolate across
        interp_mask = fx > drop_band(1) & fx < drop_band(2);

        % Ensure at least some points are in the interpolation region
        if any(interp_mask)
            % Use linear interpolation across the dropout band
            fx_clean = fx(~interp_mask);
            psd_clean = psd(~interp_mask);
            psd_interp = interp1(fx_clean, psd_clean, fx, 'linear', 'extrap');
            psd(interp_mask) = psd_interp(interp_mask);
        end

       
        idx = fx > f_range(1) & fx < f_range(2);

        fx = fx(idx);
        psd = psd(idx);

        % Take the logarithm of frequency and PSD
        log_f = log10(fx);
        log_psd = log10(psd);

        % Fit a straight line to the log-log plot using polyfit
        p = polyfit(log_f, log_psd, 1);

        % The slope of the line corresponds to the spectral slope
        spectral_slope(i) = p(1);
        
    end
end
