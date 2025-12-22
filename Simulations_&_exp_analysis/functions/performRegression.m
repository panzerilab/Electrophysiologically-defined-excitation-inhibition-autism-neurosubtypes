function [B, R2, y_pred, y_fit, bint, R2_CI] = performRegression(X, y, n_boot)

    % Add interaction term if exactly two predictors
    if size(X, 2) == 2
        X = [X, X(:,1) .* X(:,2)];
    end

    % Add intercept
    X_aug = [ones(size(X,1),1), X];

    % Perform OLS regression
    [B, bint, ~, ~, stats] = regress(y, X_aug);

    % Predictions
    y_pred = X_aug * B;

    % R² from stats
    R2 = stats(1);

    % Fit line for visualization
    p = polyfit(y, y_pred, 1);
    y_fit = polyval(p, y);

    % --- Bootstrap for R² confidence interval ---
    n_boot = n_boot;
    n = size(X_aug, 1);
    R2_boot = zeros(n_boot, 1);

    for i = 1:n_boot
        idx = randsample(n, n, true); % bootstrap sample
        Xb = X_aug(idx, :);
        yb = y(idx);
        yb_pred = Xb * (Xb \ yb);
        SS_res = sum((yb - yb_pred).^2);
        SS_tot = sum((yb - mean(yb)).^2);
        R2_boot(i) = 1 - (SS_res / SS_tot);
    end

    R2_CI = prctile(R2_boot, [2.5 97.5]); % 95% CI
    %R2_CI = std(R2_boot);

end
