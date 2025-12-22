function [spectral_slope, freq, amplitude, bw, num_peaks] = compute_spectral_slope_fooof(signals, fs, f_range)

    num_signals = size(signals, 2);
    spectral_slope = zeros(1, num_signals);
    freq = zeros(1, num_signals);
    amplitude = zeros(1, num_signals);
    bw = zeros(1, num_signals);
    num_peaks = zeros(1, num_signals); 

    params.tapers = [3, 5];
    settings.peak_width_limits = [2, 70];
    %settings.peak_width_limits = [2, 8];
    settings.max_n_peaks = 2;
    settings.aperiodic_mode = 'fixed';
    settings.peak_threshold = 1.;
    settings.verbose = false;

    for i = 1:num_signals
        clear fooof_results
        signal = signals(:, i);

        % Compute power spectrum
        [S, f] = mtspectrumc(signal, params);
        psd = S';
        fx = f * fs;

        % Interpolate to avoid drop around 55 Hz (European powerline)
        drop_band = [46, 54];  % Hz range to interpolate across
        interp_mask = fx > drop_band(1) & fx < drop_band(2);

        % Ensure at least some points are in the interpolation region
        if any(interp_mask)
            % Use linear interpolation across the dropout band
            fx_clean = fx(~interp_mask);
            psd_clean = psd(~interp_mask);
            psd_interp = interp1(fx_clean, psd_clean, fx, 'linear', 'extrap');
            psd(interp_mask) = psd_interp(interp_mask);
        end


        % Fit FOOOF
        fooof_results = fooof(fx, psd, f_range, settings);
        spectral_slope(i) = fooof_results.aperiodic_params(end);

        % Store number of peaks
        if ~isempty(fooof_results.peak_params)
            num_peaks(i) = size(fooof_results.peak_params, 1);

            % Find the peak with maximum amplitude
            [~, max_idx] = max(fooof_results.peak_params(:, 2));
            peak = fooof_results.peak_params(max_idx, :);  % [freq, amp, bw]

            if peak(1) > 5
                freq(i) = peak(1);
                amplitude(i) = peak(2);
                bw(i) = peak(3);
            else
                freq(i) = 0;
                amplitude(i) = 0;
                bw(i) = 0;
            end
        else
            num_peaks(i) = 0;
            freq(i) = 0;
            amplitude(i) = 0;
            bw(i) = 0;
        end
    end
end
