function basicFxn_drawCorr_permutation( ...
    X_allIter_allSubj, Y_allIter_allSubj, colors, x_ticks, y_ticks, ...
    x_ticklabels, y_ticklabels, flag_zeroMean, flag_plotIdvdCI, ...
    flag_plotUnikSymbol, str_title, markers_allSubj, nIter)

% =========================================================================
% basicFxn_drawCorr_permutation (REVISED)
%
% - NO type_corr / type_tail inputs
% - Always reports BOTH Pearson (r) and Spearman (rho), two-tailed
% - Partial (control loc) reported via:
%     (A) demean within loc then corr
%     (B) partialcorr with location dummy coding
% - Per-location correlations also reported (r and rho)
% - Permutation p-values (two-tailed) computed on MEDIANS
% - Bootstrap CIs (row bootstrap over subjects)
% - Linear regression fits (point + bootstrap CI bands) per location and global
% =========================================================================

%% ---------------- settings ----------------
wd_border = 3;
sz_ticks  = 40;
sz_marker = 18;

CI_plot_point = 0.68;   % per-point CI bars (X/Y across iterations) if flag_plotIdvdCI==1
CI_level      = 0.95;   % CIs reported for correlations and regression bands
nPerm         = 1e4;
nBoot         = 5000;

nSamplesLine  = 200;    % resolution for regression line plotting

%% ---------------- validate / reshape ----------------
assert(ndims(X_allIter_allSubj) == 3 && ndims(Y_allIter_allSubj) == 3, ...
    'ALERT: X_allIter_allSubj and Y_allIter_allSubj must be 3D arrays.');

nCond = size(colors, 1);              % treat "condition" as "location"
nSubj = numel(markers_allSubj);

X_allIter_allSubj = fxn_reshape(X_allIter_allSubj, nSubj, nCond, 'X_allIter_allSubj');
Y_allIter_allSubj = fxn_reshape(Y_allIter_allSubj, nSubj, nCond, 'Y_allIter_allSubj');

[nIter2, nSubj2, nCond2] = size(X_allIter_allSubj);
assert(all(size(Y_allIter_allSubj) == [nIter2, nSubj2, nCond2]), ...
    'ALERT: X and Y must match after reshape.');
assert(nIter2 == nIter, 'ALERT: Input nIter=%d but inferred nIter=%d.', nIter, nIter2);
assert(nSubj2 == nSubj && nCond2 == nCond, 'ALERT: Reshaped X has wrong nSubj/nCond.');

%% ---------------- subject-level medians and per-point CI (from iterations) ----------------
% Outputs of getCI here are [nSubj x nCond]
[X_med, X_lb, X_ub] = getCI(X_allIter_allSubj, 1, 1, CI_plot_point);
[Y_med, Y_lb, Y_ub] = getCI(Y_allIter_allSubj, 1, 1, CI_plot_point);

%% ---------------- optional centering ----------------
switch flag_zeroMean
    case 1 % subtract mean across observers within each location
        shiftX = mean(X_med, 1, 'omitnan');    % 1 x nCond
        shiftY = mean(Y_med, 1, 'omitnan');    % 1 x nCond

        X_obs = X_med - shiftX;
        Y_obs = Y_med - shiftY;

        X_lb_plot = X_lb - shiftX;
        X_ub_plot = X_ub - shiftX;
        Y_lb_plot = Y_lb - shiftY;
        Y_ub_plot = Y_ub - shiftY;

    case 2 % subtract mean across locations within each observer
        shiftX = mean(X_med, 2, 'omitnan');    % nSubj x 1
        shiftY = mean(Y_med, 2, 'omitnan');    % nSubj x 1

        X_obs = X_med - shiftX;
        Y_obs = Y_med - shiftY;

        X_lb_plot = X_lb - shiftX;
        X_ub_plot = X_ub - shiftX;
        Y_lb_plot = Y_lb - shiftY;
        Y_ub_plot = Y_ub - shiftY;

    otherwise % 0
        X_obs = X_med;
        Y_obs = Y_med;

        X_lb_plot = X_lb;  X_ub_plot = X_ub;
        Y_lb_plot = Y_lb;  Y_ub_plot = Y_ub;
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

        % Plot CI for each dot
        if flag_plotIdvdCI == 1
            plot([X_lb_plot(iSubj,iCond), X_ub_plot(iSubj,iCond)], [y0 y0], '-', ...
                'Color', ones(1,3)*0.5, 'LineWidth', 1.5, 'HandleVisibility','off');
            plot([x0 x0], [Y_lb_plot(iSubj,iCond), Y_ub_plot(iSubj,iCond)], '-', ...
                'Color', ones(1,3)*0.5, 'LineWidth', 1.5, 'HandleVisibility','off');
        end

        % Plot each dot
        plot(x0, y0, mk, ...
            'MarkerFaceColor', faceColor, ...
            'MarkerEdgeColor', edgeColor, ...
            'MarkerSize', sz_marker, ...
            'LineWidth', wd_border);
    end % iSubj
end % iCond

%% ---------------- Correlation analyses (Pearson + Spearman, two-tailed) ----------------
pct = 100*[(1-CI_level)/2, 1-(1-CI_level)/2];

xvec = X_obs(:);
yvec = Y_obs(:);

% location dummies for partialcorr
Z = repmat(1:nCond, nSubj, 1);
D = dummyvar(Z(:));
D = D(:, 1:end-1);

%% (1A) Partial via demean-within-loc + corr (OBS)
X_res = X_obs - mean(X_obs, 1, 'omitnan');
Y_res = Y_obs - mean(Y_obs, 1, 'omitnan');

r_partial_demean_obs   = safeCorr(X_res(:), Y_res(:), 'Pearson');
rho_partial_demean_obs = safeCorr(X_res(:), Y_res(:), 'Spearman');

%% (1B) Partial via partialcorr + dummies (OBS)
r_partial_pc_obs   = partialcorr(xvec, yvec, D, 'type', 'Pearson',  'rows', 'complete');
rho_partial_pc_obs = partialcorr(xvec, yvec, D, 'type', 'Spearman', 'rows', 'complete');

%% (2) Per-location correlations (OBS)
r_loc_obs   = nan(1,nCond);
rho_loc_obs = nan(1,nCond);
for iCond = 1:nCond
    r_loc_obs(iCond)   = safeCorr(X_obs(:,iCond), Y_obs(:,iCond), 'Pearson');
    rho_loc_obs(iCond) = safeCorr(X_obs(:,iCond), Y_obs(:,iCond), 'Spearman');
end

%% ---------------- [NHST] Permutation tests (two-tailed) ----------------
% (A) Partial via demean: permute subject labels of Y as blocks
r_partial_demean_perm   = nan(nPerm,1);
rho_partial_demean_perm = nan(nPerm,1);

parfor iPerm = 1:nPerm
    idx = randperm(nSubj);
    Yp = Y_obs(idx,:);
    Yp_res = Yp - mean(Yp, 1, 'omitnan');

    r_partial_demean_perm(iPerm)   = safeCorr(X_res(:), Yp_res(:), 'Pearson');
    rho_partial_demean_perm(iPerm) = safeCorr(X_res(:), Yp_res(:), 'Spearman');
end

pperm_r_partial_demean   = fxn_perm_pval(r_partial_demean_perm,   r_partial_demean_obs,   'both');
pperm_rho_partial_demean = fxn_perm_pval(rho_partial_demean_perm, rho_partial_demean_obs, 'both');

% (B) Partial via partialcorr: permute Y by subject blocks
r_partial_pc_perm   = nan(nPerm,1);
rho_partial_pc_perm = nan(nPerm,1);

parfor iPerm = 1:nPerm
    idx = randperm(nSubj);
    Yp = Y_obs(idx,:);
    yp = Yp(:);

    r_partial_pc_perm(iPerm)   = partialcorr(xvec, yp, D, 'type', 'Pearson',  'rows', 'complete');
    rho_partial_pc_perm(iPerm) = partialcorr(xvec, yp, D, 'type', 'Spearman', 'rows', 'complete');
end

pperm_r_partial_pc   = fxn_perm_pval(r_partial_pc_perm,   r_partial_pc_obs,   'both');
pperm_rho_partial_pc = fxn_perm_pval(rho_partial_pc_perm, rho_partial_pc_obs, 'both');

% (C) Per-location: permute subject labels of Y as blocks
r_loc_perm   = nan(nPerm,nCond);
rho_loc_perm = nan(nPerm,nCond);

parfor iPerm = 1:nPerm
    idx = randperm(nSubj);
    Yp = Y_obs(idx,:);
    for iCond = 1:nCond
        r_loc_perm(iPerm,iCond)   = safeCorr(X_obs(:,iCond), Yp(:,iCond), 'Pearson');
        rho_loc_perm(iPerm,iCond) = safeCorr(X_obs(:,iCond), Yp(:,iCond), 'Spearman');
    end
end

pperm_r_loc   = nan(1,nCond);
pperm_rho_loc = nan(1,nCond);
for iCond = 1:nCond
    pperm_r_loc(iCond)   = fxn_perm_pval(r_loc_perm(:,iCond),   r_loc_obs(iCond),   'both');
    pperm_rho_loc(iCond) = fxn_perm_pval(rho_loc_perm(:,iCond), rho_loc_obs(iCond), 'both');
end

%% ---------------- Regression x-grids (compute once) ----------------
% Global x-grid
xminG = min(X_obs(:), [], 'omitnan');
xmaxG = max(X_obs(:), [], 'omitnan');
if xminG == xmaxG
    x_global = linspace(xminG-1e-6, xmaxG+1e-6, nSamplesLine);
else
    x_global = linspace(xminG, xmaxG, nSamplesLine);
end

% Per-location x-grids
x_loc = nan(nCond, nSamplesLine);
for iCond = 1:nCond
    xminL = min(X_obs(:,iCond), [], 'omitnan');
    xmaxL = max(X_obs(:,iCond), [], 'omitnan');
    if ~isfinite(xminL) || ~isfinite(xmaxL), continue; end
    if xminL == xmaxL
        x_loc(iCond,:) = linspace(xminL-1e-6, xmaxL+1e-6, nSamplesLine);
    else
        x_loc(iCond,:) = linspace(xminL, xmaxL, nSamplesLine);
    end
end

%% ---------------- Point-estimate regression fits (on medians) ----------------
[~, yfit_global] = local_linfit_yhat(X_obs(:), Y_obs(:), x_global);

yfit_loc = nan(nCond, nSamplesLine);
for iCond = 1:nCond
    if any(isnan(x_loc(iCond,:))), continue; end
    [~, yfit_loc(iCond,:)] = local_linfit_yhat(X_obs(:,iCond), Y_obs(:,iCond), x_loc(iCond,:));
end

%% ---------------- [Bootstrap] ONE loop for everything (coherent resamples) ----------------
r_partial_demean_boot   = nan(nBoot,1);
rho_partial_demean_boot = nan(nBoot,1);

r_partial_pc_boot       = nan(nBoot,1);
rho_partial_pc_boot     = nan(nBoot,1);

r_loc_boot              = nan(nBoot,nCond);
rho_loc_boot            = nan(nBoot,nCond);

yfit_global_boot        = nan(nBoot, nSamplesLine);
yfit_loc_boot           = nan(nBoot, nSamplesLine, nCond);

parfor iBoot = 1:nBoot
    idx = randi(nSubj, [1 nSubj]);
    Xb = X_obs(idx, :);
    Yb = Y_obs(idx, :);

    % (A) partial demean + corr
    Xb_res = Xb - mean(Xb, 1, 'omitnan');
    Yb_res = Yb - mean(Yb, 1, 'omitnan');
    r_partial_demean_boot(iBoot)   = safeCorr(Xb_res(:), Yb_res(:), 'Pearson');
    rho_partial_demean_boot(iBoot) = safeCorr(Xb_res(:), Yb_res(:), 'Spearman');

    % (B) partialcorr + dummies
    xb = Xb(:);
    yb = Yb(:);

    Zb = repmat(1:nCond, nSubj, 1);
    Db = dummyvar(Zb(:));
    Db = Db(:, 1:end-1);

    r_partial_pc_boot(iBoot)   = partialcorr(xb, yb, Db, 'type', 'Pearson',  'rows', 'complete');
    rho_partial_pc_boot(iBoot) = partialcorr(xb, yb, Db, 'type', 'Spearman', 'rows', 'complete');

    % (C) per-location corr
    for iCond = 1:nCond
        r_loc_boot(iBoot,iCond)   = safeCorr(Xb(:,iCond), Yb(:,iCond), 'Pearson');
        rho_loc_boot(iBoot,iCond) = safeCorr(Xb(:,iCond), Yb(:,iCond), 'Spearman');
    end

    % (D) regression fits: store yhat curves
    [~, yhatG] = local_linfit_yhat(xb, yb, x_global);
    yfit_global_boot(iBoot,:) = yhatG;

    for iCond = 1:nCond
        if any(isnan(x_loc(iCond,:))), continue; end
        [~, yhatL] = local_linfit_yhat(Xb(:,iCond), Yb(:,iCond), x_loc(iCond,:));
        yfit_loc_boot(iBoot,:,iCond) = yhatL;
    end
end

%% ---------------- Bootstrap CIs (corr) ----------------
CI_r_partial_demean   = prctile(r_partial_demean_boot(~isnan(r_partial_demean_boot)), pct);
CI_rho_partial_demean = prctile(rho_partial_demean_boot(~isnan(rho_partial_demean_boot)), pct);

CI_r_partial_pc       = prctile(r_partial_pc_boot(~isnan(r_partial_pc_boot)), pct);
CI_rho_partial_pc     = prctile(rho_partial_pc_boot(~isnan(rho_partial_pc_boot)), pct);

CI_r_loc   = nan(nCond,2);
CI_rho_loc = nan(nCond,2);
for iCond = 1:nCond
    tmp = r_loc_boot(:,iCond);  tmp = tmp(~isnan(tmp));
    if ~isempty(tmp), CI_r_loc(iCond,:) = prctile(tmp, pct); else, CI_r_loc(iCond,:) = [NaN NaN]; end

    tmp = rho_loc_boot(:,iCond); tmp = tmp(~isnan(tmp));
    if ~isempty(tmp), CI_rho_loc(iCond,:) = prctile(tmp, pct); else, CI_rho_loc(iCond,:) = [NaN NaN]; end
end

%% ---------------- Bootstrap CIs (regression bands) ----------------
yfit_global_lb = prctile(yfit_global_boot, pct(1), 1);
yfit_global_ub = prctile(yfit_global_boot, pct(2), 1);

yfit_loc_lb = nan(nCond, nSamplesLine);
yfit_loc_ub = nan(nCond, nSamplesLine);
for iCond = 1:nCond
    if all(isnan(yfit_loc_boot(:,:,iCond)), 'all'), continue; end
    yfit_loc_lb(iCond,:) = prctile(yfit_loc_boot(:,:,iCond), pct(1), 1);
    yfit_loc_ub(iCond,:) = prctile(yfit_loc_boot(:,:,iCond), pct(2), 1);
end

%% ---------------- Plot regression bands + lines ----------------
% Global band + line
patch([x_global fliplr(x_global)], [yfit_global_lb fliplr(yfit_global_ub)], ...
    ones(1,3)*0.6, 'EdgeColor','none', 'FaceAlpha', 0.20, 'HandleVisibility','off');
plot(x_global, yfit_global, '-k', 'LineWidth', wd_border*1.6, 'HandleVisibility','off');

% Per-location band + line
for iCond = 1:nCond
    if any(isnan(x_loc(iCond,:))) || any(isnan(yfit_loc(iCond,:)))
        continue;
    end

    % solid if permutation p(rho) < .1, else dashed
    if ~isnan(pperm_rho_loc(iCond)) && pperm_rho_loc(iCond) < 0.1
        ls = '-';
    else
        ls = '--';
    end

    patch([x_loc(iCond,:) fliplr(x_loc(iCond,:))], ...
        [yfit_loc_lb(iCond,:) fliplr(yfit_loc_ub(iCond,:))], ...
        colors(iCond,:), 'EdgeColor','none', 'FaceAlpha', 0.12, 'HandleVisibility','off');

    % plot(x_loc(iCond,:), yfit_loc(iCond,:), ls, 'Color', colors(iCond,:), ...
    %     'LineWidth', wd_border, 'HandleVisibility','off');
end

%% ---------------- strings ----------------
str_cross0_r = '*';
str_cross0_rho = '*';
if CI_r_partial_demean(1)*CI_r_partial_demean(2) < 0, str_cross0_r = ''; end
if CI_rho_partial_demean(1)*CI_rho_partial_demean(2) < 0, str_cross0_rho = ''; end

str_partial_demean = sprintf( ...
    'Partial (demean+corr):  r%s=%.2f [%.2f, %.2f], p=%.3f, \\rho%s=%.2f [%.2f, %.2f], p=%.3f' , ...
    str_cross0_r, r_partial_demean_obs, CI_r_partial_demean(1), CI_r_partial_demean(2), pperm_r_partial_demean, ...
    str_cross0_rho, rho_partial_demean_obs, CI_rho_partial_demean(1), CI_rho_partial_demean(2), pperm_rho_partial_demean);

str_cross0_r = '*';
str_cross0_rho = '*';
if CI_r_partial_pc(1)*CI_r_partial_pc(2) < 0, str_cross0_r = ''; end
if CI_rho_partial_pc(1)*CI_rho_partial_pc(2) < 0, str_cross0_rho = ''; end

str_partial_pc = sprintf( ...
    'Partial (partialcorr):  r%s=%.2f [%.2f, %.2f], p=%.3f, \\rho%s=%.2f [%.2f, %.2f], p=%.3f' , ...
    str_cross0_r, r_partial_pc_obs, CI_r_partial_pc(1), CI_r_partial_pc(2), pperm_r_partial_pc, ...
    str_cross0_rho, rho_partial_pc_obs, CI_rho_partial_pc(1), CI_rho_partial_pc(2), pperm_rho_partial_pc);

str_perLoc = sprintf('Corr per location:\n');
for iCond = 1:nCond
    % if mod(iCond-1,2) == 0
    %     str_perLoc = [str_perLoc sprintf('\n')];
    % else
    %     str_perLoc = [str_perLoc sprintf('   ')];
    % end

    str_cross0_r = '*';
    str_cross0_rho = '*';
    if CI_r_loc(iCond,1)*CI_r_loc(iCond,2) < 0, str_cross0_r = ''; end
    if CI_rho_loc(iCond,1)*CI_rho_loc(iCond,2) < 0, str_cross0_rho = ''; end

    str_perLoc = [str_perLoc sprintf('Cond#%d: r%s=%.2f [%.2f, %.2f] (p=%.3f), \\rho%s=%.2f [%.2f, %.2f] (p=%.3f)\n', ...
        iCond, ...
        str_cross0_r,  r_loc_obs(iCond),   CI_r_loc(iCond,1),   CI_r_loc(iCond,2),   pperm_r_loc(iCond), ...
        str_cross0_rho, rho_loc_obs(iCond), CI_rho_loc(iCond,1), CI_rho_loc(iCond,2), pperm_rho_loc(iCond))];
end %iCond

%% ---------------- axes formatting ----------------
if ~isnan(x_ticks), xticks(x_ticks); xlim(x_ticks([1 end])); end
if ~isnan(y_ticks), yticks(y_ticks); ylim(y_ticks([1 end])); end
if ~isnan(x_ticklabels), xticklabels(x_ticklabels); end
if ~isnan(y_ticklabels), yticklabels(y_ticklabels); end

axis square
ax = gca;
ax.XAxis.FontSize = sz_ticks;
ax.YAxis.FontSize = sz_ticks;
ax.LineWidth = wd_border;

%% ---------------- print text blocks (bottom-center stacked) ----------------
text(ax, 0.5, 0.02, str_partial_demean, 'Units','normalized', ...
    'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
    'FontSize', 12, 'Interpreter','tex', 'Clipping','off');

text(ax, 0.5, 0.06, str_partial_pc, 'Units','normalized', ...
    'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
    'FontSize', 12, 'Interpreter','tex', 'Clipping','off');

% text(ax, 0.5, 0.12, str_perLoc, 'Units','normalized', ...
%     'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
%     'FontSize', 10, 'Interpreter','tex', 'Clipping','off');

title(sprintf('%s\n%s\n%s\n%s\n', str_title, str_partial_demean, str_partial_pc, str_perLoc), 'FontSize', 14);

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