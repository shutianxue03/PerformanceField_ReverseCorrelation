function [PearsonR, SpearmanRho] = basicFxn_drawCorr_permutation( ...
    X_allIter_allSubj, Y_allIter_allSubj, colors, flag_UseRUseRho, x_ticks, y_ticks, ...
    x_ticklabels, y_ticklabels, flag_zeroMean, flag_plotIdvdCI, ...
    flag_plotUnikSymbol, str_title, markers_allSubj, nIter)

% =========================================================================
% basicFxn_drawCorr_permutation
% - Always reports BOTH Pearson (r) and Spearman (rho), two-tailed
% - Partial (control loc) reported via:
% (A) demean within loc then corr
% (B) partialcorr with location dummy coding
% - Per-location correlations also reported (r and rho)
% - Permutation p-values (two-tailed) computed on MEDIANS
% - Bootstrap CIs (row bootstrap over subjects)
% - Linear regression fits (point + bootstrap CI bands) per location and global

% =========================================================================

%% ---------------- settings ----------------
wd_border = 3;
sz_ticks = 40;
sz_marker = 18;

CI68 = .68; % per-point CI bars (X/Y across iterations) if flag_plotIdvdCI==1
CI95 = .95; % CIs reported for correlations and regression bands
nPerm = 1e4;
nBoot = 1e4;
seedPerm = 1;
seedBoot = 2;
nSamplesLine = 200;
sigMark = @(lb,ub) ternary((lb*ub>0), '*', '');

%% ---------------- validate / reshape ----------------
assert(ndims(X_allIter_allSubj) == 3 && ndims(Y_allIter_allSubj) == 3, ...
    'ALERT: X_allIter_allSubj and Y_allIter_allSubj must be 3D arrays.');

nCond = size(colors, 1); % treat "condition" as "location"
nSubj = numel(markers_allSubj);

X_allIter_allSubj = fxn_reshape(X_allIter_allSubj, nSubj, nCond, 'X_allIter_allSubj');
Y_allIter_allSubj = fxn_reshape(Y_allIter_allSubj, nSubj, nCond, 'Y_allIter_allSubj');

[nIter2, nSubj2, nCond2] = size(X_allIter_allSubj);
assert(all(size(Y_allIter_allSubj) == [nIter2, nSubj2, nCond2]), ...
    'ALERT: X and Y must match after reshape.');
assert(nIter2 == nIter, 'ALERT: Input nIter=%d but inferred nIter=%d.', nIter, nIter2);
assert(nSubj2 == nSubj && nCond2 == nCond, 'ALERT: Reshaped X has wrong nSubj/nCond.');

%% ---------------- subject-level medians and per-point CI (from iterations) ----------------
[X_med, X_lb, X_ub, X_sem_neg, X_sem_pos] = getCI(X_allIter_allSubj, 1, 1, CI68);
[Y_med, Y_lb, Y_ub, Y_sem_neg, Y_sem_pos] = getCI(Y_allIter_allSubj, 1, 1, CI68);

%% ---------------- optional centering ----------------
switch flag_zeroMean
    case 1 % subtract mean across observers within each location
        shiftX = mean(X_med, 1, 'omitnan');
        shiftY = mean(Y_med, 1, 'omitnan');

        X_obs = X_med - shiftX;
        Y_obs = Y_med - shiftY;

        X_lb_plot = X_lb - shiftX;
        X_ub_plot = X_ub - shiftX;
        Y_lb_plot = Y_lb - shiftY;
        Y_ub_plot = Y_ub - shiftY;

    case 2 % subtract mean across locations within each observer
        shiftX = mean(X_med, 2, 'omitnan');
        shiftY = mean(Y_med, 2, 'omitnan');

        X_obs = X_med - shiftX;
        Y_obs = Y_med - shiftY;

        X_lb_plot = X_lb - shiftX;
        X_ub_plot = X_ub - shiftX;
        Y_lb_plot = Y_lb - shiftY;
        Y_ub_plot = Y_ub - shiftY;

    otherwise % 0
        X_obs = X_med;
        Y_obs = Y_med;

        X_lb_plot = X_lb; X_ub_plot = X_ub;
        Y_lb_plot = Y_lb; Y_ub_plot = Y_ub;
end

%% ---------------- Precompute vectors + dummy controls for partialcorr ----------------
xvec = X_obs(:);
yvec = Y_obs(:);

Z = repmat(1:nCond, nSubj, 1);
D = dummyvar(Z(:));
D = D(:, 1:end-1);

% Fast bootstrap of D by indexing subject blocks
rowIx_perSubj = reshape(1:(nSubj*nCond), nSubj, nCond); % [nSubj x nCond]

%% ---------------- Observed correlations ----------------
X_res = X_obs - mean(X_obs, 1, 'omitnan');
Y_res = Y_obs - mean(Y_obs, 1, 'omitnan');

r_partial_demean_obs = safeCorr(X_res(:), Y_res(:), 'Pearson');
rho_partial_demean_obs = safeCorr(X_res(:), Y_res(:), 'Spearman');

r_partial_pc_obs = partialcorr(xvec, yvec, D, 'type', 'Pearson', 'rows', 'complete');
rho_partial_pc_obs = partialcorr(xvec, yvec, D, 'type', 'Spearman', 'rows', 'complete');

r_obs_allCond = nan(1,nCond);
rho_obs_allCond = nan(1,nCond);
for iCond = 1:nCond
    r_obs_allCond(iCond) = safeCorr(X_obs(:,iCond), Y_obs(:,iCond), 'Pearson');
    rho_obs_allCond(iCond) = safeCorr(X_obs(:,iCond), Y_obs(:,iCond), 'Spearman');
end

%% ---------------- [NHST] Permutation tests (ONE parfor) ----------------
r_partial_demean_perm = nan(nPerm,1);
rho_partial_demean_perm = nan(nPerm,1);

r_partial_pc_perm = nan(nPerm,1);
rho_partial_pc_perm = nan(nPerm,1);

r_perm_allCond = nan(nPerm,nCond);
rho_perm_allCond = nan(nPerm,nCond);

% Pregenerate subj indices
rng(seedPerm, 'twister');
indRand_allPerm = zeros(nPerm, nSubj, 'uint16');
for iPerm = 1:nPerm
    indRand_allPerm(iPerm, :) = uint16(randperm(nSubj));
end

% timePermStart = tic;
parfor iPerm = 1:nPerm
    indRandPerm = double(indRand_allPerm(iPerm, :));
    Yrand = Y_obs(indRandPerm,:);

    % (A) demean + corr null
    YrandRes = Yrand - mean(Yrand, 1, 'omitnan');
    r_partial_demean_perm(iPerm) = safeCorr(X_res(:), YrandRes(:), 'Pearson');
    rho_partial_demean_perm(iPerm) = safeCorr(X_res(:), YrandRes(:), 'Spearman');

    % (B) partialcorr + dummies null
    YrandVec = Yrand(:);
    r_partial_pc_perm(iPerm) = partialcorr(xvec, YrandVec, D, 'type', 'Pearson', 'rows', 'complete');
    rho_partial_pc_perm(iPerm) = partialcorr(xvec, YrandVec, D, 'type', 'Spearman', 'rows', 'complete');

    % (C) per-location null
    for iCond = 1:nCond
        r_perm_allCond(iPerm,iCond) = safeCorr(X_obs(:,iCond), Yrand(:,iCond), 'Pearson');
        rho_perm_allCond(iPerm,iCond) = safeCorr(X_obs(:,iCond), Yrand(:,iCond), 'Spearman');
    end
end % iPerm
% timePermEnd = toc(timePermStart);
% fprintf('Permutation tests completed in %.1f seconds (%.2f min).\n', timePermEnd, timePermEnd/60);

pperm_r_partial_demean = fxn_perm_pval(r_partial_demean_perm, r_partial_demean_obs, 'both');
pperm_rho_partial_demean = fxn_perm_pval(rho_partial_demean_perm, rho_partial_demean_obs, 'both');

pperm_r_partial_pc = fxn_perm_pval(r_partial_pc_perm, r_partial_pc_obs, 'both');
pperm_rho_partial_pc = fxn_perm_pval(rho_partial_pc_perm, rho_partial_pc_obs, 'both');

pperm_r_allCond = nan(1,nCond);
pperm_rho_allCond = nan(1,nCond);
for iCond = 1:nCond
    pperm_r_allCond(iCond) = fxn_perm_pval(r_perm_allCond(:,iCond), r_obs_allCond(iCond), 'both');
    pperm_rho_allCond(iCond) = fxn_perm_pval(rho_perm_allCond(:,iCond), rho_obs_allCond(iCond), 'both');
end

%% ---------------- Regression x-grids (compute once) ----------------
xminG = min(X_obs(:), [], 'omitnan');
xmaxG = max(X_obs(:), [], 'omitnan');
if ~isfinite(xminG) || ~isfinite(xmaxG)
    x_global = nan(1,nSamplesLine);
else
    if xminG == xmaxG
        x_global = linspace(xminG-1e-6, xmaxG+1e-6, nSamplesLine);
    else
        x_global = linspace(xminG, xmaxG, nSamplesLine);
    end
end

x_perCond = nan(nCond, nSamplesLine);
for iCond = 1:nCond
    xminL = min(X_obs(:,iCond), [], 'omitnan');
    xmaxL = max(X_obs(:,iCond), [], 'omitnan');
    if ~isfinite(xminL) || ~isfinite(xmaxL), continue; end
    if xminL == xmaxL
        x_perCond(iCond,:) = linspace(xminL-1e-6, xmaxL+1e-6, nSamplesLine);
    else
        x_perCond(iCond,:) = linspace(xminL, xmaxL, nSamplesLine);
    end
end

%% ---------------- Point-estimate regression fits (on medians) ----------------
[~, yfit_global] = local_linfit_yhat(X_obs(:), Y_obs(:), x_global);

yfit_allCond = nan(nCond, nSamplesLine);
for iCond = 1:nCond
    if any(isnan(x_perCond(iCond,:))), continue; end
    [~, yfit_allCond(iCond,:)] = local_linfit_yhat(X_obs(:,iCond), Y_obs(:,iCond), x_perCond(iCond,:));
end

%% ---------------- Bootstrapping ----------------
r_partial_demean_allBoot = nan(nBoot,1);
rho_partial_demean_allBoot = nan(nBoot,1);

r_partial_pc_allBoot = nan(nBoot,1);
rho_partial_pc_allBoot = nan(nBoot,1);

r_boot_allCond = nan(nBoot,nCond);
rho_boot_allCond = nan(nBoot,nCond);

yfit_global_allBoot = nan(nBoot, nSamplesLine);
yfit_boot_allCond = nan(nBoot, nSamplesLine, nCond);

% Pregenerate subj indices
rng(seedBoot, 'twister');
indRand_allBoot = randi(nSubj, [nBoot, nSubj], 'uint16');

% timeBootStart = tic;    
parfor iBoot = 1:nBoot
    indRandBoot = double(indRand_allBoot(iBoot, :));

    Xrand = X_obs(indRandBoot, :);
    Yrand = Y_obs(indRandBoot, :);

    % (A) partial demean + corr
    XrandRes = Xrand - mean(Xrand, 1, 'omitnan');
    YrandRes = Yrand - mean(Yrand, 1, 'omitnan');
    r_partial_demean_allBoot(iBoot) = safeCorr(XrandRes(:), YrandRes(:), 'Pearson');
    rho_partial_demean_allBoot(iBoot) = safeCorr(XrandRes(:), YrandRes(:), 'Spearman');

    % (B) partialcorr + dummies (reuse D blocks)
    XrandVec = Xrand(:);
    YrandVec = Yrand(:);

    rows = rowIx_perSubj(indRandBoot, :);
    rows = rows(:);
    Db = D(rows, :);

    r_partial_pc_allBoot(iBoot) = partialcorr(XrandVec, YrandVec, Db, 'type', 'Pearson', 'rows', 'complete');
    rho_partial_pc_allBoot(iBoot) = partialcorr(XrandVec, YrandVec, Db, 'type', 'Spearman', 'rows', 'complete');

    % (C) per-location corr
    for iCond = 1:nCond
        r_boot_allCond(iBoot,iCond) = safeCorr(Xrand(:,iCond), Yrand(:,iCond), 'Pearson');
        rho_boot_allCond(iBoot,iCond) = safeCorr(Xrand(:,iCond), Yrand(:,iCond), 'Spearman');
    end % iCond

    % (D) regression fits
    [~, yhatG] = local_linfit_yhat(XrandVec, YrandVec, x_global);
    yfit_global_allBoot(iBoot,:) = yhatG;
    for iCond = 1:nCond
        if any(isnan(x_perCond(iCond,:))), continue; end
        [~, yfit_boot_allCond(iBoot,:,iCond)] = local_linfit_yhat(Xrand(:,iCond), Yrand(:,iCond), x_perCond(iCond,:));
    end % iCond

end % iBoot
% timeBootEnd = toc(timeBootStart);
% fprintf('Bootstrapping completed in %.1f seconds (%.2f min).\n', timeBootEnd, timeBootEnd/60);

%% ---------------- Bootstrap CIs (corr) ----------------
[r_partial_demean_med, r_partial_demean_lb, r_partial_demean_ub] = getCI(r_partial_demean_allBoot, 1, 1, CI95);
[rho_partial_demean_med, rho_partial_demean_lb, rho_partial_demean_ub] = getCI(rho_partial_demean_allBoot, 1, 1, CI95);

[r_partial_pc_med, r_partial_pc_lb, r_partial_pc_ub] = getCI(r_partial_pc_allBoot, 1, 1, CI95);
[rho_partial_pc_med, rho_partial_pc_lb, rho_partial_pc_ub] = getCI(rho_partial_pc_allBoot, 1, 1, CI95);

r_allCond_lb = nan(nCond,1);
r_allCond_ub = nan(nCond,1);
rho_allCond_lb = nan(nCond,1);
rho_allCond_ub = nan(nCond,1);

for iCond = 1:nCond
    [r_allCond_med, r_allCond_lb(iCond), r_allCond_ub(iCond)] = getCI(r_boot_allCond(:,iCond), 1, 1, CI95);
    [rho_allCond_med, rho_allCond_lb(iCond), rho_allCond_ub(iCond)] = getCI(rho_boot_allCond(:,iCond), 1, 1, CI95);
end

%% ---------------- Bootstrap CIs (regression bands, 68% CI) ----------------
% Global fits
[~, yfit_global_lb, yfit_global_ub] = getCI(yfit_global_allBoot, 1, 1, CI68);

% Per-condition fits
yfit_lb_allLoc = nan(nCond, nSamplesLine);
yfit_ub_allLoc = nan(nCond, nSamplesLine);
for iCond = 1:nCond
    tmp = squeeze(yfit_boot_allCond(:,:,iCond)); % [nBoot x nSamplesLine]
    if all(isnan(tmp), 'all'), continue; end
    [~, yfit_lb_allLoc(iCond,:), yfit_ub_allLoc(iCond,:)] = getCI(tmp, 1, 1, CI68);
end

%% ---------------- PLOT: scatter of medians ----------------
figure('Position', [0 200 1000 1000]);
hold on; box on;

for iCond = 1:nCond
    for iSubj = 1:nSubj

        if flag_plotUnikSymbol
            mk = markers_allSubj{iSubj};
        else
            mk = 'o';
        end

        if ismember(mk, {'+','x','*','.'})
            faceColor = 'none';
            edgeColor = colors(iCond,:);
        else
            faceColor = colors(iCond,:);
            edgeColor = 'w';
        end

        x0 = X_obs(iSubj, iCond);
        y0 = Y_obs(iSubj, iCond);

        if flag_plotIdvdCI == 1
            % errorbar([X_lb_plot(iSubj,iCond), X_ub_plot(iSubj,iCond)], [y0 y0], '-', 'Color', ones(1,3)*0.5, 'LineWidth', 1.5, 'HandleVisibility','off');
            % plot([x0 x0], [Y_lb_plot(iSubj,iCond), Y_ub_plot(iSubj,iCond)], '-', 'Color', ones(1,3)*0.5, 'LineWidth', 1.5, 'HandleVisibility','off');
            errorbar(X_med(iSubj,iCond), Y_med(iSubj,iCond), Y_sem_neg(iSubj,iCond), Y_sem_pos(iSubj,iCond), 'vertical', '-', 'Color', ones(1,3)*0.5, 'LineWidth', 1.5, 'HandleVisibility','off', 'CapSize', 0)
        end

        plot(x0, y0, mk, ...
            'MarkerFaceColor', faceColor, ...
            'MarkerEdgeColor', edgeColor, ...
            'MarkerSize', sz_marker, ...
            'LineWidth', wd_border);
    end
end

%% ---------------- Plot regression bands + lines ----------------
if ~any(isnan(x_global)) && ~any(isnan(yfit_global))
    patch([x_global, fliplr(x_global)], [yfit_global_lb, fliplr(yfit_global_ub)], ones(1,3)*0.6, 'EdgeColor','none', 'FaceAlpha', 0.20, 'HandleVisibility','off');
    plot(x_global, yfit_global, '-k', 'LineWidth', wd_border*1.6, 'HandleVisibility','off');
end

for iCond = 1:nCond
    if any(isnan(x_perCond(iCond,:))) || any(isnan(yfit_allCond(iCond,:)))
        continue;
    end

    if ~isnan(pperm_rho_allCond(iCond)) && pperm_rho_allCond(iCond) < 0.1
        ls = '-';
    else
        ls = '--';
    end

    % Plot 68% CI
    patch([x_perCond(iCond,:) fliplr(x_perCond(iCond,:))], [yfit_lb_allLoc(iCond,:) fliplr(yfit_ub_allLoc(iCond,:))], colors(iCond,:), 'EdgeColor','none', 'FaceAlpha', 0.12, 'HandleVisibility','off');

    % plot the median fits (optional)
    % plot(x_perCond(iCond,:), yfit_allCond(iCond,:), ls, 'Color', colors(iCond,:), 'LineWidth', wd_border, 'HandleVisibility','off');
end

%% ---------------- strings ----------------
str_partial_demean = sprintf( ...
    'Partial (demean+corr): r%s=%.2f [%.2f, %.2f], p=%.3f, \\rho%s=%.2f [%.2f, %.2f], p=%.3f' , ...
    sigMark(r_partial_demean_lb, r_partial_demean_ub), ...
    r_partial_demean_med, r_partial_demean_lb, r_partial_demean_ub, pperm_r_partial_demean, ...
    sigMark(rho_partial_demean_lb, rho_partial_demean_ub), ...
    rho_partial_demean_med, rho_partial_demean_lb, rho_partial_demean_ub, pperm_rho_partial_demean);

str_partial_pc = sprintf( ...
    'Partial (partialcorr): r%s=%.2f [%.2f, %.2f], p=%.3f, \\rho%s=%.2f [%.2f, %.2f], p=%.3f' , ...
    sigMark(r_partial_pc_lb, r_partial_pc_ub), ...
    r_partial_pc_med, r_partial_pc_lb, r_partial_pc_ub, pperm_r_partial_pc, ...
    sigMark(rho_partial_pc_lb, rho_partial_pc_ub), ...
    rho_partial_pc_med, rho_partial_pc_lb, rho_partial_pc_ub, pperm_rho_partial_pc);

str_perLoc = sprintf('Corr per location:\n');
for iCond = 1:nCond
    str_perLoc = [str_perLoc sprintf( ...
        'Cond#%d: r%s=%.2f [%.2f, %.2f] (p=%.3f), \\rho%s=%.2f [%.2f, %.2f] (p=%.3f)\n', ...
        iCond, ...
        sigMark(r_allCond_lb(iCond), r_allCond_ub(iCond)), ...
        r_obs_allCond(iCond), r_allCond_lb(iCond), r_allCond_ub(iCond), pperm_r_allCond(iCond), ...
        sigMark(rho_allCond_lb(iCond), rho_allCond_ub(iCond)), ...
        rho_obs_allCond(iCond), rho_allCond_lb(iCond), rho_allCond_ub(iCond), pperm_rho_allCond(iCond))];
end

%% ---------------- axes formatting ----------------
% if ~isnan(x_ticks), xticks(x_ticks); xlim(x_ticks([1 end])); end
% if ~isnan(y_ticks), yticks(y_ticks); ylim(y_ticks([1 end])); end
% if ~isnan(x_ticklabels), xticklabels(x_ticklabels); end
% if ~isnan(y_ticklabels), yticklabels(y_ticklabels); end

axis square
ax = gca;
ax.XAxis.FontSize = sz_ticks;
ax.YAxis.FontSize = sz_ticks;
ax.LineWidth = wd_border;

%% ---------------- print stats ----------------
switch flag_UseRUseRho
    case 'useR'
        str_print = sprintf('Partial r=%.2f [%.2f, %.2f]' , r_partial_pc_med, r_partial_pc_lb, r_partial_pc_ub);
    case 'useRho'
        str_print = sprintf('Partial \\rho=%.2f [%.2f, %.2f]' , rho_partial_pc_med, rho_partial_pc_lb, rho_partial_pc_ub);
end


y_str = 0.02; % Figure 6: 0.02
text(ax, 0.5, y_str, str_print, 'Units', 'normalized', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    'FontSize', 45, 'Color', 'k', 'Interpreter', 'tex', 'Clipping', 'off');

title(sprintf('%s\n%s\n%s\n%s\n', str_title, str_partial_demean, str_partial_pc, str_perLoc), 'FontSize', 14);

%% Save CI for automatic CI range calcuation in CorrAsym
PearsonR = [r_partial_pc_med, r_partial_pc_lb, r_partial_pc_ub];
SpearmanRho = [rho_partial_pc_med, rho_partial_pc_lb, rho_partial_pc_ub];

end

%% ========================= helper: safeCorr =========================
function r = safeCorr(x,y,typeName)
x = x(:); y = y(:);
ok = ~isnan(x) & ~isnan(y);
x = x(ok); y = y(ok);
if numel(x) < 3 || numel(unique(x)) < 2 || numel(unique(y)) < 2
    r = NaN; return;
end
r = corr(x, y, 'Type', typeName, 'Rows', 'complete');
end

%% ========================= helper: permutation p-value =========================
function p = fxn_perm_pval(nullStats, obsStat, type_tail)
nullStats = nullStats(~isnan(nullStats));
n = numel(nullStats);

switch lower(type_tail)
    case {'both','two','two-sided','twosided'}
        p = (1 + sum(abs(nullStats) >= abs(obsStat))) / (n + 1);
    case {'right','greater'}
        p = (1 + sum(nullStats >= obsStat)) / (n + 1);
    case {'left','less'}
        p = (1 + sum(nullStats <= obsStat)) / (n + 1);
    otherwise
        error('Unknown type_tail: %s', type_tail);
end
end

%% ========================= helper: linear fit -> yhat on grid =========================
function [beta, yhat] = local_linfit_yhat(x, y, xgrid)
x = x(:); y = y(:);
ok = ~isnan(x) & ~isnan(y);
x = x(ok); y = y(ok);

if numel(x) < 3 || numel(unique(x)) < 2
    beta = [NaN NaN];
    yhat = nan(size(xgrid));
    return;
end

beta = polyfit(x, y, 1);
yhat = polyval(beta, xgrid);
end

%% ===== helper: simple ternary =====
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end