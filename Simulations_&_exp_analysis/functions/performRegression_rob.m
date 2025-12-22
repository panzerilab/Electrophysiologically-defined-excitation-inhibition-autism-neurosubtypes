function [B, R2, y_pred, y_fit, bint, R2_CI, r, p_value, R2_boot] = performRegression(X, y, n_boot)
    
    % Perform robust regression (robustfit adds intercept internally)
    [B, stats] = robustfit(X, y,"huber"); 
    % Use a custom tuning constant (larger = more permissive)

    % Predicted values
    X_aug = [ones(size(X,1),1), X];
    y_pred = X_aug * B;

    % Compute pseudo-R²
    SS_res = sum(stats.resid .^ 2);
    SS_tot = sum((y - mean(y)) .^ 2);
    R2 = 1 - (SS_res / SS_tot);

    [r,p_value] = corr(y, y_pred, 'Type', 'Spearman');

    % Line fit for visualization
    p = polyfit(y, y_pred, 1);
    y_fit = polyval(p, y);

    % Confidence intervals for coefficients
    bint = [B - 1.96 * stats.se, B + 1.96 * stats.se];

    % --- Bootstrap confidence interval for R² ---
    n_boot = n_boot;
    n = size(X, 1);
    R2_boot = zeros(n_boot, 1);

    for i = 1:n_boot
        idx = randsample(n, n, true);  % resample with replacement
        Xb = X(idx, :);
        yb = y(idx);

        % Robust regression on bootstrap sample
        [B_b, stats_b] = robustfit(Xb, yb,"huber");
        Xb_aug = [ones(size(Xb,1),1), Xb];
        yb_pred = Xb_aug * B_b;

        SS_res_b = sum((yb - yb_pred).^2);
        SS_tot_b = sum((yb - mean(yb)).^2);
        R2_boot(i) = 1 - (SS_res_b / SS_tot_b);

    end

    R2_CI = prctile(R2_boot, [2.5 97.5]);  % 95% CI for pseudo-R²
    %R2_CI = std(R2_boot);

end
