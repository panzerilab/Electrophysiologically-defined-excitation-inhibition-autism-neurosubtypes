function [spectral_slope, freq, amplitude, bw, num_peaks] = compute_spectral_slope_fooof(signals, fs, f_range)
    % signals: Matrix of input signals (each signal is a column)
    % fs: Sampling frequency
    % f_range: Frequency range for analysis

    num_signals = size(signals, 2);
    spectral_slope = zeros(1, num_signals);
    freq = zeros(1, num_signals);
    amplitude = zeros(1, num_signals);
    bw = zeros(1, num_signals);
    num_peaks = zeros(1, num_signals); % Initialize array for number of peaks

    params.taper = [3, 5];
    settings.peak_width_limits = [1, 70];
    %settings.peak_width_limits = [2, 8];
    settings.max_n_peaks = 2;
    settings.aperiodic_mode = 'fixed'; 
    settings.peak_threshold = 0.5;
    settings.verbose = false;

    for i = 1:num_signals
        clear fooof_results
        signal = signals(:, i);

        % Compute power spectrum
        [S, f] = mtspectrumc(signal, params);
        psd = S';
        fx = f * fs;

        % Fit FOOOF
        fooof_results = fooof(fx, psd, f_range, settings);
        spectral_slope(i) = fooof_results.aperiodic_params(end);

        % Store number of peaks
        if ~isempty(fooof_results.peak_params)
            num_peaks(i) = size(fooof_results.peak_params, 1);

            % Find the peak with maximum amplitude
            [~, max_idx] = max(fooof_results.peak_params(:, 2));
            peak = fooof_results.peak_params(max_idx, :);  % [freq, amp, bw]

            if peak(1) > 30
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
