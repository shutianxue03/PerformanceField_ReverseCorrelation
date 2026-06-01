function fxn_drawBars( ...
    data_allIter_allSubj, ref, colors, x_ticklabels, y_ticks, y_ticklabels, ...
    flag_plotIDVD, flag_plotDiff, str_title, sz_fig, nIter, markers_allSubj, sz_text, wd)

% fxn_drawBars
% - Omnibus RM one-way ANOVA (per iter -> CI; perm p via within-subj label shuffle on median)
% - Pairwise posthoc paired t-tests for ALL pairs (per iter -> CI; perm p via sign-flip on median diffs)
% - Bonferroni correction for pairwise perm p-values
% - Only print significant pairs (p_perm < pThresh_print) in the title

%% ---- settings ----
% sz_text = 22; % behav, sep: 25; tunC: 22; nLL: 35;  NOM params: 25
% wd = 3; % behav: 3; sep: 2; tunC: 3

wd_bar = .5;
sz_title = 10;
sz_marker_idvd = 10;

CI95 = .95;
CI68 = .68;          % 68% for plotting (visual readability)
nPerm = 1e4; % permutation, to derive p value
nBoot = 1e4; % bootstrapping, to derive point and interval estimate
seedPerm = 1;
seedBoot = 2;

fxn_getG = @(n, dz) dz*(1-3/(4*(n-1)-1)); % Obtain Hedges' g, to correct for small sample size

fprintf('\n\n *** %s ***\n', str_title);

%% ---- shape checks: force to [nIter x nSubj x nCond] ----
assert(ndims(data_allIter_allSubj) == 3, 'ALERT: data_allIter_allSubj must be 3D arrays.');

nCond = size(colors, 1);
nSubj = numel(markers_allSubj);

data_allIter_allSubj = fxn_reshape(data_allIter_allSubj, nSubj, nCond, 'data_allIter_allSubj');

[nIter2, nSubj2, nCond2] = size(data_allIter_allSubj);
assert(nIter2 == nIter, 'ALERT: Input nIter=%d but inferred nIter=%d from X.', nIter, nIter2);
assert(nSubj2 == nSubj && nCond2 == nCond, 'ALERT: Reshaped matrix has wrong nSubj/nCond.');

pairs  = nchoosek(1:nCond, 2);
nPairs = size(pairs, 1);

%% Obtain median and CI of the data
[data_med_allSubj, ~, ~] = getCI(data_allIter_allSubj, 1, 1, CI95); % [nSubj x nCond] median over iters

%% ----------------- [NHST] Permutation p-values for ANOVA and paired contrasts (subject-level) -----------------

% ---------------- (1) Analysis of variance ----------------
[Fvalue_obs, ~] = rm_oneway(data_med_allSubj);
df_ANOVA_num = nCond - 1;
df_ANOVA_den = (nSubj - 1) * (nCond - 1);
Fvalue_allPerm = nan(nPerm, 1);

% Pregenerate subj indices
rng(seedPerm, 'twister');
indRand_allPerm = zeros(nPerm, nSubj, nCond, 'uint16');
parfor iPerm = 1:nPerm
    for iSubj = 1:nSubj
        indRand_allPerm(iPerm, iSubj, :) = uint16(randperm(nCond));
    end
end

% Generate the null distribution of F-values by shuffling condition labels within each subject
parfor iPerm = 1:nPerm
    data_perPerm = data_med_allSubj;
    for iSubj = 1:nSubj
        indRandPerm = double(squeeze(indRand_allPerm(iPerm, iSubj, :)));
        data_perPerm(iSubj, :) = data_perPerm(iSubj, indRandPerm);  % shuffle condition labels within subject
    end
    Fvalue_allPerm(iPerm) = rm_oneway(data_perPerm); % rm_oneway can return [F,eta2p]; 1 output -> F
end % iPerm

pperm_ANOVA = (1 + sum(Fvalue_allPerm >= Fvalue_obs)) / (nPerm + 1);

% ---------------- (2) Planned pairwise contrasts ----------------
t_obs_allPairs  = nan(1, nPairs);   % observed |t| per pair
CohenD_obs_allPairs  = nan(1, nPairs);   % Effect size 1: Cohen's d = mean(diff)/sd(diff)
HedgeG_obs_allPairs  = CohenD_obs_allPairs;   % Effect size 2: Hedges' g = J*dz, J=1-3/(4*(n-1)-1)
pperm_allPairs  = nan(1, nPairs);   % permutation p-values per pair (two-tailed)

% Generate the null distribution of t-values by sign-flipping the paired differences within each subject
rng(seedPerm, 'twister');
sgnMat = int8((rand(nSubj, nPerm) > 0.5) * 2 - 1); % [nSubj x nPerm]

for iPair = 1:nPairs
    iCondA = pairs(iPair, 1);
    iCondB = pairs(iPair, 2);

    diffk = data_med_allSubj(:, iCondA) - data_med_allSubj(:, iCondB);  % [nSubj x 1]
    indOK = ~isnan(diffk);
    diffk = diffk(indOK);
    nk    = numel(diffk);

    % if nk < 2 || std(diffk, 0) == 0
    %     t_obs_allPairs(iPair) = NaN;
    %     CohenD_obs_allPairs(iPair) = NaN;
    %     HedgeG_obs_allPairs(iPair) = NaN;
    %     pperm_allPairs(iPair) = NaN;
    %     continue;
    % end

    sd_diff = std(diffk, 0);
    denom   = sd_diff / sqrt(nk);

    % observed stats
    t_obs_allPairs(iPair) = abs(mean(diffk) / denom);  % |t|
    CohenD_obs_allPairs(iPair) = mean(diffk) / sd_diff;     % dz (signed)
    HedgeG_obs_allPairs(iPair) = fxn_getG(nk, CohenD_obs_allPairs(iPair));

    % sign-flip null: only numerator changes
    sgn      = double(sgnMat(indOK, :));                 % [nk x nPerm]
    num_perm = mean(diffk .* sgn, 1);            % [1 x nPerm]
    t_perm   = abs(num_perm / denom);            % [1 x nPerm]

    % two-tailed permutation p-value (+1 correction)
    pperm_allPairs(iPair) = (1 + sum(t_perm >= t_obs_allPairs(iPair))) / (nPerm + 1);
end % iPair

% ---------------- (3) Comparison to the reference (per condition) ----------------
if ~isnan(ref)
    % Initialize outputs with the variable names you implied
    tRef_obs_allCond    = nan(1, nCond);   % observed |t| for (cond - ref)
    CohenD_Ref_obs_allCond    = nan(1, nCond);   % observed Cohen's d for (cond - ref)
    HedgeG_Ref_obs_allCond = CohenD_Ref_obs_allCond;
    pperm_ref_allCond   = nan(1, nCond);   % permutation p-values per condition (two-tailed)

    % Pre-generate sign flips once (consistent null draws)
    rng(seedPerm+1e3, 'twister');
    sgnMat_ref = int8((rand(nSubj, nPerm) > 0.5) * 2 - 1);   % [nSubj x nPerm]

    % Generate the null distribution of t-values
    for iCond = 1:nCond
        diffk = data_med_allSubj(:, iCond) - ref;      % [nSubj x 1]
        indOK = ~isnan(diffk);
        diffk = diffk(indOK);
        nk    = numel(diffk);

        % if nk < 2 || std(diffk, 0) == 0
        %     tRef_obs_allCond(iCond)  = NaN;
        %     dzRef_obs_allCond(iCond)  = NaN;
        %     gzRef_obs_allCond(iCond)  = NaN;
        %     pperm_ref_allCond(iCond) = NaN;
        %     continue;
        % end

        sd_diff = std(diffk, 0);
        denom   = sd_diff / sqrt(nk);

        % observed
        tRef_obs_allCond(iCond) = abs(mean(diffk) / denom);  % |t|
        CohenD_Ref_obs_allCond(iCond) = mean(diffk) / sd_diff;     % Cohen's d (signed)
        HedgeG_Ref_obs_allCond(iCond) = fxn_getG(nk, CohenD_Ref_obs_allCond(iCond));     % Hedges' g (signed)

        % sign-flip null
        sgn      = double(sgnMat_ref(indOK, :));                % [nk x nPerm]
        num_perm = mean(diffk .* sgn, 1);               % [1 x nPerm]
        t_perm   = abs(num_perm / denom);               % [1 x nPerm]

        pperm_ref_allCond(iCond) = (1 + sum(t_perm >= tRef_obs_allCond(iCond))) / (nPerm + 1);
    end % iCond
end

%% ---------------- [Bootstrap] CIs for planned contrasts (row bootstrap over subjects) --------------

CI_diffPair = 1 - 0.05;   % planned -> 95%
CI_refPair  = 1 - 0.05;   % planned -> 95%

% Empty placeholders
groupAve_allBoot = nan(nBoot, nCond);
eta2p_allBoot = nan(nBoot, 1);
diffPair_allBoot = nan(nBoot, nPairs); % difference among conditions
CohenDPair_allBoot = diffPair_allBoot; % Cohen's d
HedgeGPair_allBoot = CohenDPair_allBoot; % Hedges' g

diffRef_allBoot = nan(nBoot, nCond); % difference between each cond and the ref
CohenDRef_allBoot = diffRef_allBoot;
HedgeGRef_allBoot = CohenDRef_allBoot;

% Pregenerate subj indices
rng(seedBoot, 'twister');
indRand_allBoot = randi(nSubj, [nBoot, nSubj]);

parfor iBoot = 1:nBoot
    % resample rows (subjects) with replacement
    % indResampled = randi(nSubj, [1, nSubj]);
    indRandBoot = double(indRand_allBoot(iBoot, :));
    dataRand = data_med_allSubj(indRandBoot, :);   % [nSubj x nCond] bootstrap sample

    % Group average (for plotting)
    groupAve_allBoot(iBoot, :) = mean(dataRand, 1);  % bootstrap mean per condition

    % Analysis of variance
    [~, eta2p_allBoot(iBoot)] = rm_oneway(dataRand);

    % Pairwise planned differences
    for iPair = 1:nPairs
        iCondA = pairs(iPair, 1);
        iCondB = pairs(iPair, 2);
        diffk = dataRand(:, iCondA) - dataRand(:, iCondB);
        diffPair_allBoot(iBoot, iPair) = mean(diffk);
        CohenDPair_allBoot(iBoot, iPair) = mean(diffk, 'omitnan')/std(diffk);
        HedgeGPair_allBoot(iBoot, iPair) = fxn_getG(numel(diffk), CohenDPair_allBoot(iBoot, iPair));

    end % iPair

    % Comparisons to reference
    if ~isnan(ref)
        for iCond = 1:nCond
            diffk = dataRand(:, iCond) - ref;
            diffRef_allBoot(iBoot, iCond) = mean(diffk);
            CohenDRef_allBoot(iBoot, iCond) = mean(diffk)/std(diffk);
            HedgeGRef_allBoot(iBoot, iCond) = fxn_getG(numel(diffk), CohenDRef_allBoot(iBoot, iCond));
        end % iCond
    end
end % iBoot

%% Obtain point and interval estimates
% 68% CI for plotting
[groupAve_med, ~, ~, groupAve_sem_neg68, groupAve_sem_pos68] = getCI(groupAve_allBoot, 1, 1, CI68);
% 95% CI for reporting
[~, groupAve_lb95, groupAve_ub95] = getCI(groupAve_allBoot, 1, 1, CI95);

[eta2p_med, eta2p_lb, eta2p_ub] = getCI(eta2p_allBoot, 1, 1, CI95);
[diffPair_med, diffPair_lb, diffPair_ub] = getCI(diffPair_allBoot, 1, 1, CI_diffPair); % for reporting
[~, ~, ~, diffPair_sem_neg, diffPair_sem_pos] = getCI(diffPair_allBoot, 1, 1, CI68); % for plotting
% [CohenDPair_med, CohenDPair_lb, CohenDPair_ub] = getCI(CohenDPair_allBoot, 1, 1, CI_diffPair);
[HedgeGPair_med, HedgeGPair_lb, HedgeGPair_ub] = getCI(HedgeGPair_allBoot, 1, 1, CI_diffPair);
[~, diffRef_lb, diffRef_ub] = getCI(diffRef_allBoot, 1, 1, CI_refPair);
% [CohenDRef_med, CohenDRef_lb, CohenDRef_ub] = getCI(CohenDRef_allBoot, 1, 1, CI_refPair);
[HedgeGRef_med, HedgeGRef_lb, HedgeGRef_ub] = getCI(HedgeGRef_allBoot, 1, 1, CI_refPair);

%% Strings handling

% (1) Analysis of variance
str_ANOVA = sprintf('ANOVA: F(%d,%d)=%.2f, p=%.3f, eta2p=%.2f [%.2f, %.2f]', df_ANOVA_num, df_ANOVA_den, Fvalue_obs, pperm_ANOVA, eta2p_med, eta2p_lb, eta2p_ub);

% (2) Pairwise comparisons
str_diffPair = sprintf('Planned pairwise contrasts (%.1f%%CI)\n', CI_diffPair*100);
for iPair = 1:nPairs
    str_cross0 = '';
    if diffPair_lb(iPair)*diffPair_ub(iPair)>0, str_cross0 = '*'; end
    str_diffPair = [str_diffPair, sprintf('%s-%s%s=%.2f [%.2f, %.2f] | t=%.2f, p=%.3f, g=%.2f [%.2f, %.2f]\n', ...
        x_ticklabels{pairs(iPair,1)}, x_ticklabels{pairs(iPair,2)}, str_cross0, ...
        diffPair_med(iPair), diffPair_lb(iPair), diffPair_ub(iPair), ...
        t_obs_allPairs(iPair), pperm_allPairs(iPair), HedgeGPair_med(iPair), HedgeGPair_lb(iPair), HedgeGPair_ub(iPair))];
end

% (3) Compare to the reference
str_diffRef = 'No ref';
if ~isnan(ref)
    str_diffRef = sprintf('Ref=%.1f (%.2f%% CI)\n', ref, CI_refPair);

    for iCond = 1:nCond
        str_cross0 = '';
        if diffRef_lb(iCond)*diffRef_ub(iCond)>0, str_cross0 = '*'; end
        str_diffRef = [str_diffRef, sprintf('%s%s: %.2f [%.2f, %.2f] | t=%.2f, p=%.3f, g=%.2f [%.2f, %.2f]\n', ...
            x_ticklabels{iCond}, str_cross0, ...
            groupAve_med(iCond), groupAve_lb95(iCond), groupAve_ub95(iCond), ...
            tRef_obs_allCond(iCond), pperm_ref_allCond(iCond), HedgeGRef_med(iCond), HedgeGRef_lb(iCond), HedgeGRef_ub(iCond))];
    end % iCond
end

%% [Plot] Group averages
if ~isnan(sz_fig), figure('Position', [0, 200, sz_fig]); end
hold on

for iCond = 1:nCond
    % Bootstrapped group average
    bar(iCond, groupAve_med(iCond), 'FaceColor', colors(iCond,:), 'EdgeColor', colors(iCond,:), 'BarWidth', wd_bar, 'HandleVisibility', 'off');

    % Error bars for bootstrap CI
    % 95% CI
    % errorbar(iCond, groupAve_med(iCond), groupAve_sem_neg95(iCond), groupAve_sem_pos95(iCond), '.', 'color', ones(1,3)/2, 'CapSize', 0, 'linewidth', wd/2, 'HandleVisibility', 'off');
    % 68% CI
    errorbar(iCond, groupAve_med(iCond), groupAve_sem_neg68(iCond), groupAve_sem_pos68(iCond), '.', 'color', colors(iCond,:), 'CapSize', 0, 'linewidth', wd, 'HandleVisibility', 'off');

    % Fancy white overlay (purely for aesthetics)
    if groupAve_med(iCond) > 0
        errorbar(iCond, groupAve_med(iCond), groupAve_sem_neg68(iCond), 0, '.w', 'CapSize', 0, 'linewidth', wd, 'HandleVisibility', 'off');
    else
        errorbar(iCond, groupAve_med(iCond), 0, groupAve_sem_pos68(iCond), '.w', 'CapSize', 0, 'linewidth', wd, 'HandleVisibility', 'off');
    end
end % iCond

%% [Plot] Idvd medians
if flag_plotIDVD
    if nCond == 2
        xIDVD = [1.3, 1.7];
    else
        xIDVD = 1:nCond;
    end

    for iSubj = 1:nSubj
        plot(xIDVD, data_med_allSubj(iSubj,:), '-', 'color', ones(1,3)*.7, 'markerfacecolor','w', 'markeredgecolor', ones(1,3)*.7, 'markersize', sz_marker_idvd, 'linewidth', wd/1.5);
    end
end

%% [Plot] Reference
if ~isnan(ref)
    yline(ref, '--', 'color', ones(1,3)/2, 'handlevisibility', 'off', 'linewidth', wd);
end

%% [Plot] ticks, limits
xticks(1:nCond); xticklabels([])
% if ~isempty(x_ticklabels), xticklabels(x_ticklabels), end
if ~isnan(y_ticks)
    yticks(y_ticks);
    ylim(y_ticks([1, end]));
end
if ~isnan(y_ticklabels)
    yticklabels(y_ticklabels);
end

if flag_plotIDVD
    buffer = .6;
else
    buffer = .5;
end
xlim([1-buffer, nCond+buffer]);

ax = gca;
ax.XAxis.FontSize = sz_text;
ax.YAxis.FontSize = sz_text;
ax.LineWidth = wd;

%% Prepare the y-pos of the comparison line
yl = ylim;                      % [ymin ymax]
yMin = yl(1);
yMax = yl(2);

% Base height at 80% of y-axis span
yBar = yMin + 0.80 * (yMax - yMin);

% If CI would exceed yMax, nudge downward
yTop = yBar + diffPair_sem_pos;
if yTop > yMax
    yBar = yMax - diffPair_sem_pos - 0.02 * (yMax - yMin);
end

%% [Plot] "diff bar" (nCond==2 only) (must be placed after axis are set)
% Use bootstrap CI of the mean paired difference (not SEM).
if nCond == 2 && flag_plotDiff

    % --- draw horizontal line between the two conditions ---
    plot([1, 2], [yBar, yBar], 'k-', 'LineWidth', wd, 'HandleVisibility', 'off');

    % --- draw CI errorbar at the center ---
    errorbar(1.5, yBar, diffPair_sem_neg, diffPair_sem_pos, 'k-', 'LineWidth', wd, 'HandleVisibility', 'off', 'CapSize', 0);

    % Prepare x- and y-pos of the text
    xText = 1.5;                                  % left aligned near bar 1
    yPad  = 0.02 * (yMax - yMin);                  % padding in axis units
    yText = (yBar + diffPair_sem_pos) + yPad;          % above upper CI

    % If that would exceed yMax, clamp a bit
    if yText > yMax
        yText = yMax - 0.01 * (yMax - yMin);
    end

    % Prepare the string
    str_delta = sprintf('$\\Delta=%s$ [%s, %s]\n$\\mathit{p}=%.3f$', ...
        fxn_formatSignedDecimal(diffPair_med), fxn_formatSignedDecimal(diffPair_lb), fxn_formatSignedDecimal(diffPair_ub), pperm_allPairs);

    % Print
    text(xText, yText, str_delta, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', sz_text, ...
        'Color', 'k', ...
        'Interpreter', 'latex', ...
        'Clipping', 'off');
end


%% Plot p value of ANOVA results If nCond>2
if nCond > 2
    % --- draw horizontal line across conditions ---
    plot([1, nCond], [yBar, yBar], 'k-', 'LineWidth', wd, 'HandleVisibility', 'off');

    % Prepare the string with LaTeX formatting (keep everything in one math environment)
    str_ANOVA_p = sprintf('$F(%d,%d) = %.2f, \\; \\mathit{p} = %.3f$', df_ANOVA_num, df_ANOVA_den, Fvalue_obs, pperm_ANOVA);

    % Prepare x- and y-pos of the text
    xText = (1 + nCond) / 2;                                  % placed at the middle
    yPad  = 0.1 * (yMax - yMin);                  % padding in axis units
    yText = yBar + yPad;          % above upper CI

    % If that would exceed yMax, clamp a bit
    if yText > yMax
        yText = yMax - 0.01 * (yMax - yMin);
    end

    % Print
    text(xText, yText, str_ANOVA_p, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', sz_text, ...
        'Color', 'k', ...
        'Interpreter', 'latex', ...
        'Clipping', 'off');
end

%% [Plot] Title
title(sprintf('%s\n%s\n%s%s\n\n', str_title, str_ANOVA, str_diffPair, str_diffRef), 'fontsize', sz_title);


end
%% 
function str = fxn_formatSignedDecimal(val)
if isnan(val)
    str = 'NaN';
    return;
end

abs_val = abs(val);
if abs_val >= 1
    str = sprintf('%+.2f', val);
elseif abs_val >= 0.01
    str = sprintf('%+.2f', val);
else
    str = sprintf('%+.4f', val);
    str = regexprep(str, '(\.\d*?)0+$', '$1');
    str = regexprep(str, '\.$', '');
end
end
