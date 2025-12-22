function [delta, theta, alpha, beta, gamma, broad_1, broad_2] = compute_bands(signals, fs)

    % signals: Matrix of input signals (each column is a signal)
    % fs: Sampling frequency

    num_signals = size(signals, 2);
    params.tapers = [3, 5];

    % Define frequency band edges (in Hz)
    edges = [0, 4; 4, 8; 8, 12; 12, 30; 30, 70];

    % Additional broad bands
    broad_1_band = [20, 50];
    broad_2_band = [20, 80];

    % Preallocate output arrays
    delta = zeros(1, num_signals);
    theta = zeros(1, num_signals);
    alpha = zeros(1, num_signals);
    beta = zeros(1, num_signals);
    gamma = zeros(1, num_signals);
    broad_1 = zeros(1, num_signals);
    broad_2 = zeros(1, num_signals);

    for i = 1:num_signals
        signal = signals(:, i);
        
        % Compute power spectral density (PSD)
        [S, f] = mtspectrumc(signal, params);
        psd = S';
        fx = fs * f; % Frequency in Hz

        % Interpolate to fix powerline dropout (44-56 Hz)
        drop_band = [44, 56];  
        interp_mask = fx > drop_band(1) & fx < drop_band(2);
        if any(interp_mask)
            fx_clean = fx(~interp_mask);
            psd_clean = psd(~interp_mask);
            psd_interp = interp1(fx_clean, psd_clean, fx, 'linear', 'extrap');
            psd(interp_mask) = psd_interp(interp_mask);
        end

        % Sum PSD for standard bands
        for j = 1:size(edges, 1)
            idx = (fx > edges(j, 1)) & (fx <= edges(j, 2));
            band_power = sum(psd(idx));
            switch j
                case 1, delta(i) = band_power;
                case 2, theta(i) = band_power;
                case 3, alpha(i) = band_power;
                case 4, beta(i) = band_power;
                case 5, gamma(i) = band_power;
            end
        end

        % Sum PSD for broad bands
        broad_1(i) = sum(psd(fx > broad_1_band(1) & fx <= broad_1_band(2)));
        broad_2(i) = sum(psd(fx > broad_2_band(1) & fx <= broad_2_band(2)));

    end
end
