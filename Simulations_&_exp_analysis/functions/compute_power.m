function [tot_pow,gamma_pow,low_pow, freq, PSD] = compute_power(signals, fs, f_range)
    % signals: Matrix of input signals (each signal is a column)
    % fs: Sampling frequency
    % f_range: Frequency range for analysis
    num_signals = size(signals, 2);
    params.tapers=[3,5];

    for i = 1:num_signals
        
        signal = signals(:, i);
   
        [S, f] = mtspectrumc(signal,params);
        fx=fs*f;
        psd = S';
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

        freq(i,:) = fs*f;
        PSD(i,:) = psd;

        idx=fx>0 & fx<90;
        fx=fx(idx);
        psd=psd(idx);
        tot_pow(i) = sum(psd);
      
        fx=fs*f;
        psd = S';
        idx=fx>f_range(1) & fx<f_range(2);
        fx=fx(idx);
        psd=psd(idx);
        gamma_pow(i) = sum(psd);

        
        fx=fs*f;
        psd = S';
        idx= fx<f_range(1);
        fx=fx(idx);
        psd=psd(idx);
        low_pow(i) = sum(psd);

        
    end
    
end
