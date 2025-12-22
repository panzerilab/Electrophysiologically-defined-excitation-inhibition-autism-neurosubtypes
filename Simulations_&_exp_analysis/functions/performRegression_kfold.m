function [B_boot, r_boot] = performRegression_kfold(X, y, n_folds, random_seed)

    % ---- 10-fold cross-validation ----
    %n_folds = 10;
    rng(random_seed); 
    cvp = cvpartition(size(X,1), 'KFold', n_folds);
    y_pred = nan(size(y));
    r_boot = zeros(n_folds, 1);

    for f = 1:n_folds
        train_idx = training(cvp, f);
        test_idx  = test(cvp, f);

        % Robust regression on training folds
        [B_fold, ~] = robustfit(X(train_idx,:), y(train_idx), "logistic");

        % Predict on test fold
        X_test_aug = [ones(sum(test_idx),1) X(test_idx,:)];
        y_pred(test_idx) = X_test_aug * B_fold;
        
        r_boot(f) = corr(y(test_idx), X_test_aug * B_fold, 'Type','Pearson');
        B_boot(f,:) = B_fold;

    end 


end

