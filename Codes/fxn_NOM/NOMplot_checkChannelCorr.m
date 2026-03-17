%% ===== Debug: inspect design matrix conditioning =====
    % Choose the same predictor matrix used for regression
    X = reshape(e3D_tmpl_rand_norm, size(e3D_tmpl_rand_norm,1), []);
    % Remove any columns with zero variance just in case
    colSD = std(X, 0, 1);
    badCols = (colSD == 0) | isnan(colSD);
    if any(badCols)
        fprintf('Removing %d zero-variance / NaN columns out of %d total channels.\n', ...
            sum(badCols), size(X,2));
        X = X(:, ~badCols);
    end

    % Basic size
    fprintf('X size: %d trials x %d channels\n', size(X,1), size(X,2));

    % Correlation matrix across channels
    R = corr(X);

    % Singular values and condition number
    s = svd(X, 'econ');
    condX = s(1) / s(end);

    fprintf('Largest singular value: %.4g\n', s(1));
    fprintf('Smallest singular value: %.4g\n', s(end));
    fprintf('Condition number cond(X): %.4g\n', condX);

    % Also inspect X''X if you want
    XtX = X' * X;
    eigvals = eig(XtX);
    eigvals = sort(real(eigvals), 'descend');

    condXtX = eigvals(1) / eigvals(end);
    fprintf('Condition number cond(X''X): %.4g\n', condXtX);

    % Effective rank / cumulative variance of singular values
    varExplained = cumsum(s.^2) / sum(s.^2);
    n90 = find(varExplained >= 0.90, 1, 'first');
    n95 = find(varExplained >= 0.95, 1, 'first');
    fprintf('Components to explain 90%% variance: %d\n', n90);
    fprintf('Components to explain 95%% variance: %d\n', n95);

    % Plot channel correlation matrix
    figure('Position', [0 0 2e3 800]);
    subplot(2,3,1), hold on
    imagesc(R);
    axis square;
    colorbar;
    title('Channel correlation matrix');

    % Plot singular values
    subplot(2,3,2), hold on
    plot(s, 'o-');
    xlabel('Component index');
    ylabel('Singular value');
    title('Singular values of design matrix X');

    % Plot cumulative explained variance
    subplot(2,3,3), hold on
    plot(varExplained, 'o-');
    xlabel('Component index');
    ylabel('Cumulative variance explained');
    title('Cumulative variance explained by singular values');
    ylim([0 1]);

    % Optional: histogram of off-diagonal channel correlations
    mask = ~eye(size(R));
    rvals = R(mask);
    subplot(2,3,4), hold on
    histogram(rvals, 50);
    xlabel('Channel-channel correlation');
    ylabel('Count');
    title('Off-diagonal correlations across channels');

    % Average correlation of each channel with all others =====
    Xfull = reshape(e3D_tmpl_rand_norm, size(e3D_tmpl_rand_norm,1), []);
    Rfull = corr(Xfull);

    avgAbsCorr = mean(abs(Rfull - eye(size(Rfull))), 2);   % exclude self-correlation
    avgAbsCorr2D = reshape(avgAbsCorr, size(e3D_tmpl_rand_norm,2), size(e3D_tmpl_rand_norm,3));

    subplot(2,3,5), hold on
    imagesc(avgAbsCorr2D);
    axis square;
    colorbar;
    xlabel('SF channel');
    ylabel('Orientation channel');
    title('Average absolute correlation of each channel with others');

    saveas(gcf, sprintf('%s/check_ChannelCorr.jpg', nameFolder_Figures_perSubj))
    close all