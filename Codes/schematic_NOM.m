% NOM illustration figures
clear; close all; clc;

SX_RC1_setting
% Define folder for saving figures
nameFolder_Fig_NOM_Trialwise = sprintf('%s/NOM_Trialwise_%d%d', nameFolder_Figures, nORI, nSF);
if isempty(dir(nameFolder_Fig_NOM_Trialwise)), mkdir(nameFolder_Fig_NOM_Trialwise), end

nameFolder_Fig_Schematic = sprintf('%s/Schematic', nameFolder_Fig_NOM_Trialwise);
if isempty(dir(nameFolder_Fig_Schematic)), mkdir(nameFolder_Fig_Schematic), end

%% 1D probability distribution for pYES
color_shade = [0.8 0.7 0.9]; % light purple
color_cri = [0.4 0 0.6];
facealpha = .3;

for flag_stage = 1:5

    % ----------------- Stage names ---------------------------------------
    switch flag_stage
        case 1, str_stage = 'justIV';
        case 2, str_stage = 'moreIV';
        case 3, str_stage = 'oneDist';
        case 4, str_stage = 'oneDist_withCriterion';
        case 5, str_stage = 'moreDist_withCriterion';
    end

    % ----------------- Basic parameters ----------------------------------
    % Use the same IVs across stages for consistency
    rng(1); % for reproducibility
    IV_all = [1, .1, 2]; % 3 example IVs
    SDadd = 1; % additive noise
    criterion = IV_all(1)+.5; % criterion for right tail

    % Decide which IV(s) are used in this stage
    switch flag_stage
        case {1, 3, 4} % single-IV stages
            idx_IV = 1; % use the first IV (0)
        case {2, 5} % multi-IV stages
            idx_IV = 1:3; % use all three
    end

    IV = IV_all(idx_IV);
    Nmul = 1;
    sigma = sqrt((IV .* Nmul).^2 + SDadd.^2); % one sigma per IV

    % Flags for what to draw
    drawDist = ismember(flag_stage, [3 4 5]);
    drawCriterion = ismember(flag_stage, [4 5]);
    drawBandwidth = ismember(flag_stage, [3 4]); % bandwidth arrows for stages 3 & 4

    % ----------------- Set up figure -------------------------------------
    figure('Position',[200 200 800 500]); hold on;

    % For storing axis ranges
    x_min = min(IV);
    x_max = max(IV);
    y_max_global = 0; % will update when we draw PDFs

    % ----------------- Draw distributions (if needed) --------------------
    if drawDist
        nIV = numel(IV);

        % Colors for multiple distributions (subtle variations)
        baseColor = [0 0 0]; % black for curves

        for iTrial = 1:nIV
            mu_i = IV(iTrial);
            sigma_i = sigma(iTrial);

            % x-range for this Gaussian
            x_i = linspace(mu_i - 4*sigma_i, mu_i + 4*sigma_i, 1000);
            y_i = 1/(sigma_i*sqrt(2*pi)) * exp(-(x_i - mu_i).^2 / (2*sigma_i^2));

            % Plot Gaussian
            plot(x_i, y_i, 'Color', baseColor, 'LineWidth', 2);

            % Update global ranges
            x_min = min(x_min, min(x_i));
            x_max = max(x_max, max(x_i));
            y_max_global = max(y_max_global, max(y_i));

            % Vertical line at mean
            plot([mu_i mu_i], [0 max(y_i)], '-', 'color', [.5, .5, .5], 'LineWidth', 1.5);

            % Bandwidth (FWHM) as horizontal arrow (only for stages 3 & 4)
            if drawBandwidth
                ymax_i = max(y_i);
                y_half = ymax_i/2;
                FWHM = 2*sqrt(2*log(2)) * sigma_i;
                x1 = mu_i - FWHM/2;
                x2 = mu_i + FWHM/2;

                h = quiver(x1, y_half, x2 - x1, 0, 0, 'k', 'LineWidth', 1.5);
                h.MaxHeadSize = 0.03; % slim arrow head
            end

            % Shaded tail to the right of criterion (for stages 4 & 5)
            if drawCriterion
                idx_tail = x_i >= criterion;
                if any(idx_tail)
                    x_fill = [x_i(idx_tail), fliplr(x_i(idx_tail))];
                    y_fill = [zeros(1, sum(idx_tail)), fliplr(y_i(idx_tail))];

                    patch(x_fill, y_fill, color_shade, 'FaceAlpha', facealpha, 'EdgeColor', 'none');
                end
            end
        end
    end

    % ----------------- Criterion line (if needed) ------------------------
    if drawCriterion
        % use global y max (or a default if no PDFs were drawn for some reason)
        if y_max_global == 0, y_max_global = 1; end
        plot([criterion criterion], [0 .5], '--', 'Color',color_cri, 'LineWidth', 2);
    end

    % ----------------- Draw IV dots at y=0 -----------------
    plot(IV, zeros(size(IV)), 'ok', 'MarkerSize', 16, 'MarkerFaceColor', 'w', 'LineWidth', 2);

    % ----------------- Axes formatting -----------------------------------
    % If no distributions were drawn, set a reasonable y-range
    if y_max_global == 0
        y_max_global = 1;
    end

    % Add some margins on x
    x_margin = 0.5;
    x_min = -4;
    x_max = 6;
    y_max_global = .4;
    xlim([x_min - x_margin, x_max + x_margin]);
    ylim([0, y_max_global*1.1]);

    ax = gca;
    ax.Box = 'off'; % removes top & right box edges
    ax.YColor = 'none'; % hide y-axis line & ticks
    ax.XTick = []; % no x-ticks
    ax.LineWidth = 2;
    set(gca, 'FontSize', 14);
    set(findall(gcf, '-property', 'linewidth'), 'linewidth',3.5)

    % ----------------- Save figure ---------------------------------------
    saveas(gcf, sprintf('%s/NOM_stage%d_%s.jpg', nameFolder_Fig_Schematic, flag_stage, str_stage));

end % flag_stage
close all

%% 2D probability distribution for pA (not assuming independence)
V = 0.2; % mean of decision variable for both passes
sigma = 1; % SD for both passes

rhoDV1 = 0.1; % correlation for first layer
rhoDV2 = 0.6; % correlation for second layer (change as you like)

criterion = V + 0.3;
nContours = 3;

nSigma = 2.5;
% ----- Common grid -------------------------------------------------------
nGrid = 201;
x1_min = V - nSigma*sigma;
x1_max = V + nSigma*sigma;
x2_min = V - nSigma*sigma;
x2_max = V + nSigma*sigma;

[x1_grid, x2_grid] = meshgrid( ...
    linspace(x1_min, x1_max, nGrid), ...
    linspace(x2_min, x2_max, nGrid));

X = [x1_grid(:), x2_grid(:)];

% ----- First Gaussian (rhoDV1) -------------------------------------------
mu = [V, V];
Sigma1 = [ sigma^2, rhoDV1*sigma*sigma; ...
    rhoDV1*sigma*sigma, sigma^2 ];

Z1 = mvnpdf(X, mu, Sigma1);
Z1 = reshape(Z1, size(x1_grid));

% ----- Second Gaussian (rhoDV2) ------------------------------------------
Sigma2 = [ sigma^2, rhoDV2*sigma*sigma; ...
    rhoDV2*sigma*sigma, sigma^2 ];

Z2 = mvnpdf(X, mu, Sigma2);
Z2 = reshape(Z2, size(x1_grid));

% ----- Figure ------------------------------------------------------------
figure('Position',[200 200 700 600]); hold on;

% 1) Contour lines of the first joint density
[~, hC1] = contour(x1_grid, x2_grid, Z1, nContours, 'k', 'LineWidth', 1.5);
% 2) Contour lines of the second joint density (different style)
[~, hC2] = contour(x1_grid, x2_grid, Z2, nContours, '--', 'LineWidth', 1.5, 'Color', [0.3 0.3 0.3]);

% 3) Shade agreement regions (D1 > crit & D2 > crit) and (D1 < crit & D2 < crit)
% Top-right (YES–YES)
patch([criterion x1_max x1_max criterion], ...
    [criterion criterion x2_max x2_max], ...
    color_shade, 'FaceAlpha', facealpha, 'EdgeColor', 'none');

% Bottom-left (NO–NO)
patch([x1_min criterion criterion x1_min], ...
    [x2_min x2_min criterion criterion], ...
    color_shade, 'FaceAlpha', facealpha, 'EdgeColor', 'none');

% 4) Criterion lines
plot([criterion criterion], [x2_min x2_max], '--', 'Color', color_cri, 'LineWidth', 3);
plot([x1_min x1_max], [criterion criterion], '--', 'Color', color_cri, 'LineWidth', 3);

% 5) IVs
plot([V, V], [x1_min, x1_max], '-', 'color', [.5,.5, .5], 'linewidth', 2);
plot([x2_min, x2_max], [V, V], '-', 'color', [.5,.5, .5], 'linewidth', 2);
plot(V, x1_min, 'o', 'MarkerEdgeColor', [.5,.5, .5], 'MarkerFaceColor', 'w', 'MarkerSize', 15, 'linewidth', 2)
plot(x2_min, V, 'o', 'MarkerEdgeColor', [.5,.5, .5], 'MarkerFaceColor', 'w', 'MarkerSize', 15, 'linewidth', 2)

% ----- Compute pA for each rho using mvncdf ---------------
lb_YY = [criterion, criterion];
ub_YY = [inf, inf];
lb_NN = [-inf, -inf];
ub_NN = [criterion, criterion];

pA1 = mvncdf(lb_YY, ub_YY, mu, Sigma1) + mvncdf(lb_NN, ub_NN, mu, Sigma1);
pA2 = mvncdf(lb_YY, ub_YY, mu, Sigma2) + mvncdf(lb_NN, ub_NN, mu, Sigma2);

% ----- Formatting --------------------------------------------------------
xlabel('Internal variable of pass 1 ($V_{i,1}$)', 'Interpreter', 'latex');
ylabel('Internal variable of pass 2 ($V_{i,2}$)', 'Interpreter', 'latex');
axis equal;
xlim([x1_min x1_max]);
ylim([x2_min x2_max]);

ax = gca;
% ax.Box = 'off'; % removes top & right box edges
% ax.YColor = 'none'; % hide y-axis line & ticks
ax.XTick = []; % no x-ticks
ax.YTick = []; % no x-ticks
ax.LineWidth = 2;
set(gca, 'FontSize', 14);
% set(findall(gcf, '-property', 'linewidth'), 'linewidth', 3)

% title('Joint 2D Gaussian with criterion and agreement regions');

legend([hC1, hC2], ...
    {sprintf('$\\rho = %.2f, p_A=%.1f$', rhoDV1, pA1), sprintf('$\\rho = %.2f, p_A=%.1f$', rhoDV2, pA2)}, ...
    'Location', 'northwest', 'Interpreter', 'latex');

set(findall(gcf, '-property', 'fontsize'), 'fontsize', 25)
saveas(gcf, sprintf('%s/NOM_stage6_2D.jpg', nameFolder_Fig_Schematic));
close all

%% Illustration of misprediction when excluding noise parameters
close all; clc;

% ---------------- Ground-truth parameters ----------------
Nmul_true = .6; % induced (multiplicative) noise
SDadd_true = 15; % baseline (additive) noise
rho_true = 0.5; % correlation between passes

% IV range (internal variable mean)
IV = linspace(0, 150, 100); % row vector
medianIV = median(IV);

% Criterion in z-units
c_zscore = -.3;

% ---------------- Helper: full model sigma(IV) ----------------
sigma_true = sqrt((Nmul_true .* IV).^2 + SDadd_true.^2);
sigma_true = max(sigma_true, 1e-6); % avoid zero
criterion_true = c_zscore .* sigma_true + medianIV; % trial-wise criterion in IV units

% ---------------- 1D: pYES for full and reduced models ----------------
% Full model: V ~ N(IV, sigma_true^2), YES if V > criterion_true
pYES_true = 1 - normcdf(criterion_true, IV, sigma_true);

% Model without induced noise (Nmul = 0)
Nmul_noInduced = 0;
sigma_noInduced = sqrt((Nmul_noInduced .* IV).^2 + SDadd_true.^2);
sigma_noInduced = max(sigma_noInduced, 1e-6);
criterion_noInduced = c_zscore .* sigma_noInduced + medianIV;
pYES_noInduced = 1 - normcdf(criterion_noInduced, IV, sigma_noInduced);

% Model without constant noise (SDadd = 0)
SDadd_noConst = 0;
sigma_noConst = sqrt((Nmul_true .* IV).^2 + SDadd_noConst.^2);
sigma_noConst = max(sigma_noConst, 1e-6);
criterion_noConst = c_zscore .* sigma_noConst + medianIV;
pYES_noConst = 1 - normcdf(criterion_noConst, IV, sigma_noConst);

% "Model without correlation" for pYES:
% Note: correlation does NOT affect the marginal, so this is identical to pYES_true
pYES_noRho = pYES_true;

% ---------------- 2D: pA (probability of agreement) ----------------
% Agreement: both YES or both NO for two passes

pA_true = zeros(size(IV));
pA_noRho = zeros(size(IV));
pA_noInduced = zeros(size(IV));
pA_noConst = zeros(size(IV));

for iTrial = 1:numel(IV)
    mu_i = [IV(iTrial), IV(iTrial)];

    % Full model sigma and criterion at this IV
    s_true = sigma_true(iTrial);
    criterion_true_i = criterion_true(iTrial);
    Sigma_true = [s_true^2, rho_true*s_true^2; ...
        rho_true*s_true^2, s_true^2 ];

    % Sigma & criterion without induced noise (Nmul = 0)
    s_noInduced = sigma_noInduced(iTrial);
    criterion_noInd_i = criterion_noInduced(iTrial);
    Sigma_noInduced = [s_noInduced^2, rho_true*s_noInduced^2; ...
        rho_true*s_noInduced^2, s_noInduced^2 ];

    % Sigma & criterion without constant noise (SDadd = 0)
    s_noConst = sigma_noConst(iTrial);
    criterion_noConst_i = criterion_noConst(iTrial);
    Sigma_noConst = [s_noConst^2, rho_true*s_noConst^2; ...
        rho_true*s_noConst^2, s_noConst^2 ];

    % Sigma without correlation (rho = 0), same criterion as full model
    Sigma_noRho = [s_true^2, 0; ...
        0, s_true^2];
    criterion_noRho_i = criterion_true_i;

    % ----- Agreement regions for each model -----
    % Full model
    lb_YY_true = [criterion_true_i, criterion_true_i];
    ub_YY_true = [ inf, inf];
    lb_NN_true = [-inf, -inf];
    ub_NN_true = [criterion_true_i, criterion_true_i];

    P_YY = mvncdf(lb_YY_true, ub_YY_true, mu_i, Sigma_true);
    P_NN = mvncdf(lb_NN_true, ub_NN_true, mu_i, Sigma_true);
    pA_true(iTrial) = P_YY + P_NN;

    % No rho
    lb_YY_0 = [criterion_noRho_i, criterion_noRho_i];
    ub_YY_0 = [ inf, inf];
    lb_NN_0 = [-inf, -inf];
    ub_NN_0 = [criterion_noRho_i, criterion_noRho_i];

    P_YY_0 = mvncdf(lb_YY_0, ub_YY_0, mu_i, Sigma_noRho);
    P_NN_0 = mvncdf(lb_NN_0, ub_NN_0, mu_i, Sigma_noRho);
    pA_noRho(iTrial) = P_YY_0 + P_NN_0;

    % No induced noise
    lb_YY_noInd = [criterion_noInd_i, criterion_noInd_i];
    ub_YY_noInd = [ inf, inf];
    lb_NN_noInd = [-inf, -inf];
    ub_NN_noInd = [criterion_noInd_i, criterion_noInd_i];

    P_YY_noInd = mvncdf(lb_YY_noInd, ub_YY_noInd, mu_i, Sigma_noInduced);
    P_NN_noInd = mvncdf(lb_NN_noInd, ub_NN_noInd, mu_i, Sigma_noInduced);
    pA_noInduced(iTrial) = P_YY_noInd + P_NN_noInd;

    % No constant noise
    lb_YY_noConst = [criterion_noConst_i, criterion_noConst_i];
    ub_YY_noConst = [ inf, inf];
    lb_NN_noConst = [-inf, -inf];
    ub_NN_noConst = [criterion_noConst_i, criterion_noConst_i];

    P_YY_noConst = mvncdf(lb_YY_noConst, ub_YY_noConst, mu_i, Sigma_noConst);
    P_NN_noConst = mvncdf(lb_NN_noConst, ub_NN_noConst, mu_i, Sigma_noConst);
    pA_noConst(iTrial) = P_YY_noConst + P_NN_noConst;
end

% ---------------- Plotting ----------------
sz_font = 35;
figure('Position', [200 200 500 400]);hold on;
% ---- Panel 1: pYES misfit due to excluding induced/constant noise ----
% subplot(1,2,1); hold on;

plot(IV, pYES_true, 'k-', 'LineWidth', 4);
plot(IV, pYES_noInduced, 'k--', 'LineWidth', 4);
% plot(IV, pYES_noConst, 'b-', 'LineWidth', 2);

yline(0.5, '--', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off', 'LineWidth', 3);
text(IV(65), 0.42, 'pYES=0.5', 'color', [.5, .5, .5])

xlabel('Internal variable (IV)');
ylabel('Detection rate (p_{YES})');
% title('Effect of excluding additive / multiplicative noise');
% legend({'Full model', 'No induced noise', 'No constant noise'}, 'Location', 'southeast');

set(gca, 'FontSize', 12, 'LineWidth', 1.5, 'Box', 'off');
set(findall(gcf, '-property', 'fontsize'), 'fontsize', sz_font)
ylim([0 1]);
ax=gca;
ax.XTick = []; % no x-ticks
ax.YTick = []; % no x-ticks
% sgtitle('Mispredictions when excluding noise parameters');
saveas(gcf, sprintf('%s/NOM_excludeParams_pYES.jpg', nameFolder_Fig_Schematic));

% ---- Panel 2: pA misfit (here only showing full vs no correlation) ----
figure('Position', [200 200 500 400]);hold on;
% subplot(1,2,2); hold on;

plot(IV, pA_true, 'k-', 'LineWidth', 4);
plot(IV, pA_noRho, 'k--', 'LineWidth', 4);
% If you want to show these too, uncomment:
% plot(IV, pA_noInduced, 'r-', 'LineWidth', 2);
% plot(IV, pA_noConst, 'b-', 'LineWidth', 2);

yline(0.5, '--', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off', 'LineWidth', 3);
text(IV(65), 0.42, 'pA=0.5', 'color', [.5, .5, .5])
xlabel('Internal variable (IV)');
ylabel('Resp. consistency (p_A)');
% title('Effect of excluding noise correlation');

% legend({'Full model', 'No correlation'}, 'Location', 'southeast');
ax=gca;
ax.XTick = []; % no x-ticks
ax.YTick = []; % no x-ticks
set(gca, 'FontSize', 12, 'LineWidth', 1.5, 'Box', 'off');
set(findall(gcf, '-property', 'fontsize'), 'fontsize', sz_font)
ylim([0 1]);

% sgtitle('Mispredictions when excluding noise parameters');
saveas(gcf, sprintf('%s/NOM_excludeParams_pA.jpg', nameFolder_Fig_Schematic));

close all

%% Sketches for neural variability

close all; clc

% Time axis
nT = 500;
t = linspace(-pi, 4*pi, nT);   % one cycle

% ---------- 1. Neural variability (different amplitudes) ----------
y1 = sin(t);          % neuron 1
y2 = y1/2;    % neuron 2, lower amplitude

figure('Position',[100 100 350 350]), hold on
plot(t, y1, 'k-', 'LineWidth', 2);
plot(t, y2, 'k--', 'LineWidth', 2);
xlabel('Time')
ylabel('Activity')
ylim([-3, 3])
% title('Neuronal variability')
% axis tight
box off
set(gca, 'YTick', [], 'XTick', [])   % keep it schematic
axis off
saveas(gcf, sprintf('%s/NeuralVar.jpg', nameFolder_Fig_Schematic))
close all

%% Sketch for neural correlation
offset = 1;              % vertical shift
z1 = sin(t);             % neuron 1
z2 = sin(t) - offset;    % neuron 2, same fluctuations + offset


figure('Position',[100 100 350 350]), hold on
plot(t, z1, 'k-', 'LineWidth', 2);
plot(t, z2, 'k--', 'LineWidth', 2);
xlabel('Time')
ylabel('Activity')
% title('Neural correlation')
% axis tight
ylim([-4, 4])
box off
set(gca, 'YTick', [], 'XTick', [])
axis off

saveas(gcf, sprintf('%s/NeuralCorr.jpg', nameFolder_Fig_Schematic))
close all

%% Sketch for orientation and SF tuning function
iFamily_ORI = 1; paramsORI =[1 40 0];
% iFamily_ORI = 8; paramsORI =[1 40 0];
iFamily_SF = 2; paramsSF =[2 1 .3 0];
sz_wd = 10;
sz_hw = 50;
sz_hl = sz_hw;

for iFeature = 1:nFeatures
    switch iFeature
        case 1
            iFamily = iFamily_ORI;
            x = linspace(-90, 90, 1e2);
            params = paramsORI;
        case 2
            iFamily=iFamily_SF;
            x = 2.^linspace(0,2, 1e2);
            params = paramsSF;
    end
    % ---- Compute kernel ----
    y = predSFkernel(x, iFamily, params, 0);

    % ---- Create figure ----
    for flag_plotMode=1:2
        switch flag_plotMode

            case 1, str_mode = 'PeakAmp';
            case 2, str_mode = 'Bandw';
        end
        fig = figure('Position',[200+iFeature*300 200 450 200]); hold on
        plot(x, y, 'k-', 'LineWidth', sz_wd);

        axis tight
        axis off

        % ============================================================
        % Compute peak and bandwidth
        % ============================================================
        y_base = min(y);
        [y_peak, idx_peak] = max(y);
        x_peak = x(idx_peak);

        y_half = y_base + (y_peak - y_base)/2;

        idx_half = find(y >= y_half);
        x_left  = x(idx_half(1));
        x_right = x(idx_half(end));

        % draw half-height line
        % plot([x_left x_right],[y_half y_half],'--','LineWidth', sz_wd);

        % ============================================================
        % Coordinate transform for annotation arrows
        % ============================================================
        ax = gca;
        axpos = ax.Position;
        xlim_d = xlim; ylim_d = ylim;

        data2norm = @(xx,yy)[ ...
            axpos(1) + (xx - xlim_d(1))/diff(xlim_d)*axpos(3), ...
            axpos(2) + (yy - ylim_d(1))/diff(ylim_d)*axpos(4) ];

        % ============================================================
        % Add vertical arrow (peak amplitude)
        % ============================================================
        p0 = data2norm(x_peak, y_base);
        p1 = data2norm(x_peak, y_peak);
        if flag_plotMode==1
            annotation('arrow',[p0(1) p1(1)], [p0(2) p1(2)], 'Color','k','LineWidth',sz_wd, 'HeadWidth', sz_hw, 'HeadLength', sz_hl);
        end
        % ============================================================
        % Add horizontal double-arrow (bandwidth)
        % ============================================================
        b0 = data2norm(x_left,  y_half);
        b1 = data2norm(x_right, y_half);
        if flag_plotMode==2
            % Left → Right
            annotation('arrow', [b0(1) b1(1)], [b0(2) b1(2)], ...
                'Color','k','LineWidth', sz_wd, 'HeadWidth', sz_hw, 'HeadLength', sz_hl);

            % Right → Left (overlays shaft, gives you a double-headed arrow)
            annotation('arrow', [b1(1) b0(1)], [b1(2) b0(2)], ...
                'Color','k','LineWidth', sz_wd, 'HeadWidth', sz_hw, 'HeadLength', sz_hl);
        end
        % ============================================================
        % (Optional) Add vertical double-arrow spanning full y-axis
        % ============================================================
        y_top = ylim_d(2);
        y_bot = ylim_d(1);
        y_mid = mean(ylim_d);

        p_mid = data2norm(x_peak, y_mid);
        p_top = data2norm(x_peak, y_top);
        p_bot = data2norm(x_peak, y_bot);

        if flag_plotMode==1
            % upward arrow
            annotation('arrow', [p_mid(1) p_top(1)], [p_mid(2) p_top(2)], 'Color','k','LineWidth',sz_wd, 'HeadWidth', sz_hw, 'HeadLength', sz_hl);
            % downward arrow
            annotation('arrow', [p_mid(1) p_bot(1)], [p_mid(2) p_bot(2)], 'Color','k','LineWidth',sz_wd, 'HeadWidth', sz_hw, 'HeadLength', sz_hl);
        end
        % ============================================================
        % Save figure
        % ============================================================
        saveas(gcf, sprintf('%s/TuningFxn_%s_%s.jpg', nameFolder_Fig_Schematic, namesFeature{iFeature}, str_mode))
    end % flag_plotMode
end % iFeature
close all