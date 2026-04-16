%% ============================================================
% Script: simPlot_part1_criterionRecovery.m
%
% Part 1. Perfect observer (true template + IN = 0)
%
% Purpose:
% Test whether criterion can be recovered, and whether recovery is robust
% to true criterion (Cz_true) and signal contrast.
%
% Created by Shutian Xue on 04/06/2026
%% ============================================================

clear; clc; close all;

%--------------%
SX_RC1_setting;
%--------------%

%% ---------------- settings ----------------
str_part = 'Part1';

% For Part 1, true template should be used
iModelA_fit = 2;  % 2 = true / ideal template
iModelB_fit = 1;  % full model; does not matter 

nameFolder_Data = sprintf('%s/Data_%s', nameFolder_server, str_part);
nameFolder_Data_OOD = sprintf('%s/Data_OOD_%d%d', nameFolder_Data, nORI, nSF);
nameFolder_Data_NOM_Trialwise = sprintf('%s/Data_NOM_Trialwise_%d%d', nameFolder_Data, nORI, nSF);

nameFolder_Figures_part = fullfile(nameFolder_Figures, sprintf('IO_%s', str_part));
if ~exist(nameFolder_Figures_part, 'dir')
    mkdir(nameFolder_Figures_part);
end

%% ---------------- find IO folders ----------------
nameDir = dir(nameFolder_Data_OOD);
nameDir = nameDir([nameDir.isdir]);
nameDir = nameDir(~ismember({nameDir.name}, {'.', '..'}));

nFiles = numel(nameDir);

%% ---------------- storage ----------------
nameIO_all = cell(nFiles, 1);
Cz_true_all = nan(nFiles, 1);
signalCST_all = nan(nFiles, 1);

criterion_est_mean_all = nan(nFiles, 1);
criterion_est_lb_all   = nan(nFiles, 1);
criterion_est_ub_all   = nan(nFiles, 1);

criterion_DV_true_all = nan(nFiles, 1);

keepFile = false(nFiles, 1);

%% ---------------- loop over folders ----------------
for iFile = 1:nFiles

    nameIO = nameDir(iFile).name;
    nameFolder_OOD_load = fullfile(nameFolder_Data_OOD, nameIO);
    nameFolder_NOM_load = fullfile(nameFolder_Data_NOM_Trialwise, nameIO);

    fprintf('%s: %d/%d  %s\n', datetime('now'), iFile, nFiles, nameIO);

    % ----- parse signal contrast and Cz from nameIO -----
    tok = regexp(nameIO, 'IO_cN([-\d\.]+)_cG([-\d\.]+).*_Cz([-\d\.]+)_', 'tokens', 'once');
    if isempty(tok)
        fprintf('Could not parse signalCST / Cz_true from %s\n', nameIO);
        continue;
    end

    noiseCST = str2double(tok{1}) / 100;
    signalCST = str2double(tok{2}) / 100;
    Cz_true = str2double(tok{3});

    % ----- load truth -----
    nameFile_truth = fullfile(nameFolder_OOD_load, 'truth.mat');
    if ~exist(nameFile_truth, 'file')
        fprintf('Missing truth.mat: %s\n', nameFile_truth);
        continue;
    end
    truth = load(nameFile_truth);

    % ----- load fitNOM -----
    str_loadfitNOM = sprintf('n*_A%dB%d.mat', iModelA_fit, iModelB_fit);
    nameDir_fitNOM = dir(fullfile(nameFolder_NOM_load, str_loadfitNOM));
    nameDir_fitNOM = nameDir_fitNOM(~contains({nameDir_fitNOM.name}, 'min'));

    if isempty(nameDir_fitNOM)
        fprintf('No fitNOM file found in %s\n', nameFolder_NOM_load);
        continue;
    end

    nameFile_fitNOM = fullfile(nameDir_fitNOM(1).folder, nameDir_fitNOM(1).name);
    data_fitNOM = load(nameFile_fitNOM);

    % ----- get estimated criterion across iterations -----
    % adjust the variable name here if needed
    if isfield(data_fitNOM, 'params_est_allIter')
        params_est_allIter = data_fitNOM.params_est_allIter;

        % assume criterion is the last parameter in each iteration
        nIter = numel(params_est_allIter);
        criterion_est_allIter = params_est_allIter(:, end);


    elseif isfield(data_fitNOM, 'criterion_est_allIter')
        criterion_est_allIter = data_fitNOM.criterion_est_allIter(:);

    else
        fprintf('No criterion estimate found in %s\n', nameFile_fitNOM);
        continue;
    end

    criterion_est_allIter = criterion_est_allIter(isfinite(criterion_est_allIter));
    if isempty(criterion_est_allIter)
        fprintf('criterion_est_allIter empty in %s\n', nameFile_fitNOM);
        continue;
    end

    % ----- summarize across iterations -----
    [criterion_est_mean, criterion_est_lb, criterion_est_ub] = getCI(criterion_est_allIter);

    % ----- store -----
    keepFile(iFile) = true;

    nameIO_all{iFile} = nameIO;
    Cz_true_all(iFile) = Cz_true;
    signalCST_all(iFile) = signalCST;

    criterion_DV_true_all(iFile) = truth.criterion_DV_true;
    criterion_est_mean_all(iFile) = criterion_est_mean;
    criterion_est_lb_all(iFile) = criterion_est_lb;
    criterion_est_ub_all(iFile) = criterion_est_ub;
end % iFile

%% ---------------- keep valid rows only ----------------
nameIO_all = nameIO_all(keepFile);
Cz_true_all = Cz_true_all(keepFile);
signalCST_all = signalCST_all(keepFile);

criterion_DV_true_all = criterion_DV_true_all(keepFile);
criterion_est_mean_all = criterion_est_mean_all(keepFile);
criterion_est_lb_all = criterion_est_lb_all(keepFile);
criterion_est_ub_all = criterion_est_ub_all(keepFile);

nKeep = numel(Cz_true_all);
fprintf('\nKept %d valid folders.\n', nKeep);

%% ---------------- plotting setup ----------------
Cz_unique = unique(Cz_true_all);
signal_unique = unique(signalCST_all);

% color by true Cz
cmap_cz = parula(numel(Cz_unique));

% grayscale edge by signal contrast
signal_min = min(signal_unique);
signal_max = max(signal_unique);

%% ---------------- figure ----------------
figure('Position', [100 100 750 650]); hold on;

for iPt = 1:nKeep

    % face color from true Cz
    iCz = find(Cz_unique == Cz_true_all(iPt), 1, 'first');
    color_face = cmap_cz(iCz, :);

    % edge gray from signal contrast
    if signal_max == signal_min
        gray_edge = 0.4;
    else
        frac = (signalCST_all(iPt) - signal_min) / (signal_max - signal_min);
        gray_edge = 0.85 - 0.65 * frac;  % higher contrast -> darker edge
    end
    color_edge = gray_edge * [1 1 1];

    % 68% CI errorbar
    line([criterion_DV_true_all(iPt), criterion_DV_true_all(iPt)], ...
         [criterion_est_lb_all(iPt), criterion_est_ub_all(iPt)], ...
         'Color', color_edge, 'LineWidth', 1.5);

    % point
    scatter(criterion_DV_true_all(iPt), criterion_est_mean_all(iPt), 85, ...
        'MarkerFaceColor', color_face, ...
        'MarkerEdgeColor', color_edge, ...
        'LineWidth', 1.5);
end

% unity line
xy_all = [criterion_DV_true_all; criterion_est_mean_all; criterion_est_lb_all; criterion_est_ub_all];
xy_min = min(xy_all);
xy_max = max(xy_all);
plot([xy_min, xy_max], [xy_min, xy_max], 'k--', 'LineWidth', 1.5);

xlabel('True criterion (DV unit)');
ylabel('Estimated criterion (DV unit)');
title('Part 1. Perfect observer: criterion recovery');

axis square;
xlim([xy_min, xy_max]);
ylim([xy_min, xy_max]);
box on;
set(gca, 'FontSize', 12);

%% ---------------- legends ----------------
% legend for Cz color
for iCz = 1:numel(Cz_unique)
    hCz(iCz) = scatter(nan, nan, 85, ...
        'MarkerFaceColor', cmap_cz(iCz, :), ...
        'MarkerEdgeColor', 'k', ...
        'LineWidth', 1.2);
end
leg1 = legend(hCz, arrayfun(@(x) sprintf('Cz = %.1f', x), Cz_unique, 'UniformOutput', false), ...
    'Location', 'northwest');
title(leg1, 'True criterion');

% legend for signal contrast edge gray
hold on;
for iS = 1:numel(signal_unique)
    if signal_max == signal_min
        gray_edge = 0.4;
    else
        frac = (signal_unique(iS) - signal_min) / (signal_max - signal_min);
        gray_edge = 0.85 - 0.65 * frac;
    end
    color_edge = gray_edge * [1 1 1];

    hSig(iS) = scatter(nan, nan, 85, ...
        'MarkerFaceColor', [1 1 1], ...
        'MarkerEdgeColor', color_edge, ...
        'LineWidth', 1.5);
end
leg2 = legend(hSig, arrayfun(@(x) sprintf('signal = %.2f', x), signal_unique, 'UniformOutput', false), ...
    'Location', 'southeast');
title(leg2, 'Signal contrast');

% first legend: Cz color
leg1 = legend(hCz, arrayfun(@(x) sprintf('Cz = %.1f', x), Cz_unique, 'UniformOutput', false), ...
    'Location', 'northwest');
title(leg1, 'True criterion');

% second legend: signal contrast edge gray
ax1 = gca;
ax2 = axes('Position', ax1.Position, 'Color', 'none', 'Visible', 'off');
leg2 = legend(ax2, hSig, arrayfun(@(x) sprintf('signal = %.2f', x), signal_unique, 'UniformOutput', false), ...
    'Location', 'southeast');
title(leg2, 'Signal contrast');

%% ---------------- save ----------------
saveas(gcf, fullfile(nameFolder_Figures_part, 'Part1_criterionRecovery.png'));