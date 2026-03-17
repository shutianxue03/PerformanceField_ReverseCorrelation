
% Build X and y
X = reshape(e3D_tmpl_rand_norm, size(e3D_tmpl_rand_norm,1), []);
y = resp_tmpl_rand_sel(:);

% nORI = size(e3D_tmpl_rand_norm, 2);
% nSF  = size(e3D_tmpl_rand_norm, 3);

% Remove zero-variance / NaN columns
colSD = std(X, 0, 1);
goodCols = ~(colSD == 0 | isnan(colSD));
X_good = X(:, goodCols);

% Center X
muX = mean(X_good, 1);
Xc = X_good - muX;

% PCA via SVD
[U, S, V] = svd(Xc, 'econ');
scores = U * S;
singVals = diag(S);
varExplained = cumsum(singVals.^2) / sum(singVals.^2);

% Inspect cumulative variance
figure; hold on
plot(varExplained(varExplained<.95), 'o-');
xlabel(sprintf('PC index (%d# reached 95% VarExplained)', sum(varExplained<.95)));
ylabel('Cumulative variance explained');
title('PCA cumulative variance explained');
ylim([0 1]);

% Try several k values
kList = [10:5:35];
rList = nan(size(kList));
kernel_all = cell(size(kList));

for ii = 1:numel(kList)
    k = kList(ii);

    Z = scores(:, 1:k);

    % Linear regression on binary response (debugging)
    % beta_pc = [ones(size(Z,1),1), Z] \ y;
    % bZ = beta_pc(2:end);

    % Logistic regression on binary response
    [b_all, ~, ~] = glmfit(Z, y, 'binomial', 'link', 'logit');

    % intercept = b_all(1);
    bZ = b_all(2:end);   % coefficients in PC space

    % Project back to channel space
    b_channel_good = V(:,1:k) * bZ;

    b_channel_full = nan(1, numel(goodCols));
    b_channel_full(goodCols) = b_channel_good;

    kernel_pca = reshape(b_channel_full, nORI, nSF);
    kernel_all{ii} = kernel_pca;

    rList(ii) = corr(template_true(:), kernel_pca(:), 'rows', 'complete');

    fprintf('k = %d, corr(template_true, kernel_pca) = %.4f\n', k, rList(ii));

    figure;
    subplot(1,3,1);
    imagesc(filtersSF_all, filtersOri_all - 90, template_true);
    set(gca, 'YDir', 'normal');
    axis square; colorbar;
    title('True template');

    subplot(1,3,2);
    imagesc(filtersSF_all, filtersOri_all - 90, kernel_pca);
    set(gca, 'YDir', 'normal');
    axis square; colorbar;
    title(sprintf('PCA recovered, k=%d', k));

    subplot(1,3,3);
    imagesc(filtersSF_all, filtersOri_all - 90, kernel_pca - template_true);
    set(gca, 'YDir', 'normal');
    axis square; colorbar;
    title('Recovered - True');
end % ii

figure;
plot(kList, rList, 'o-');
xlabel('Number of PCs kept');
ylabel('Correlation with true template');
title('Recovery quality vs number of PCs');