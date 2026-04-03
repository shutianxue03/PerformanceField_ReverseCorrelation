%% ============================================================
% Compare derived vs true templates across criterion x lambda
% Shutian Xue  in Mar, 2026
% ============================================================

clear; clc; close all;
set(0, 'DefaultFigureVisible', 'on')

str_metric = 'template'; str_params = '';
% str_metric = 'parameters'; iParam = 3; str_params = 'Nshared_true';
% str_metric = 'parameters'; iParam = 4; str_params = 'criterion_DV_true';

flag_plotPerCond=0; % time-consuming!!

% ---------- USER SETTINGS ----------
%--------------%
SX_RC1_setting;
%--------------%
nameFolder_Data = sprintf('%s/Data_Part1', nameFolder_server) ;
nameFolder_Data_OOD = sprintf('%s/Data_OOD_%d%d', nameFolder_Data, nORI, nSF);  % Folder to save data
nameFolder_Data_NOM_Trialwise = sprintf('%s/Data_NOM_Trialwise_%d%d', nameFolder_Data, nORI, nSF);  % Folder to save data

% compIV file pattern
for iModelA_sim = [2,1];
    fprintf('\n\n ============= ModelA=%d =============\n\n', iModelA_sim)
    str_loadCompIV = sprintf('n*_A%d_compIV.mat', iModelA_sim);   % use A1 if derived template = RC-derived template
    str_loadfitNOM = sprintf('n*_A%dB1.mat', iModelA_sim);   % use A1 if derived template = RC-derived template

    % Extract file names
    nameDir = dir(nameFolder_Data_OOD);
    nameDir = nameDir([nameDir.isdir]);
    nameDir = nameDir(~ismember({nameDir.name}, {'.', '..'}));

    % Create empty folders
    nFiles = numel(nameDir);

    nameIO_allFiles = cell(nFiles, 1);
    Cz_allFiles = nan(nFiles, 1);
    lambda_allFiles = Cz_allFiles;
    rmse_allFiles = nan(nFiles, 1e3); % 1e3 is just temporary
    R2_allFiles = rmse_allFiles;
    valTrue_allFiles = Cz_allFiles;
    valData_allFiles = rmse_allFiles;
    criterion_est_allFiles = rmse_allFiles;
    Nshared_est_allFiles = rmse_allFiles;

    for iFile = 1:nFiles

        nameIO = nameDir(iFile).name;
        nameFile_full = fullfile(nameFolder_Data_NOM_Trialwise, nameIO);

        % Parse Cz and lambda from IO folder name
        tok = regexp(nameIO, 'Cz([-\d\.]+)_cont[-\d\.]+_whiten([-\d\.]+)', 'tokens', 'once');
        Cz_true = str2double(tok{1});
        lambda_true = str2double(tok{2});

        % Load the truth
        nameFile_truth = fullfile(nameFolder_Data_OOD, nameIO, 'truth.mat');
        if ~exist(nameFile_truth, 'file')
            fprintf('Missing truth.mat: %s\n', nameFile_truth);
            continue
        end
        truth = load(nameFile_truth);

        % Load compIV data
        nameDir_compDV = dir(fullfile(nameFile_full, str_loadCompIV));
        if isempty(nameDir_compDV)
            fprintf('No compIV file found in %s\n', nameFile_full);
            continue
        end
        nameFile_compIV = fullfile(nameFile_full, nameDir_compDV(1).name);
        data_compIV = load(nameFile_compIV);

        % Load fitNOM data
        nameDir_fitNOM = dir(fullfile(nameFile_full, str_loadfitNOM));
        if isempty(nameDir_fitNOM)
            fprintf('No fitNOM file found in %s\n', nameFile_full);
            continue
        end
        nameFile_fitNOM = fullfile(nameFile_full, nameDir_fitNOM(1).name);
        data_fitNOM = load(nameFile_fitNOM);

        % Load truth for nIter
        load(nameFile_truth, 'nIter')

        %% Plot all figures per step

        if flag_plotPerCond
            % Load behav
            nameFile_behav = fullfile(nameFolder_Data_OOD, nameIO, 'behavMeas.mat');
            if ~exist(nameFile_behav, 'file')
                fprintf('Missing behav.mat: %s\n', nameFile_behav);
                continue
            end
            load(nameFile_behav);

            % Load truth
            load(nameFile_truth)

            % Load compIV
            load(nameFile_compIV)

            % Load fitNOM
            load(nameFile_fitNOM)

            iLocComb=1;
            subjName = nameIO;
            iModelA = iModelA_sim;

            % Print progress

            nameFolder_Figures_perSubj = sprintf('%s/IO/%s', nameFolder_Figures, subjName);

            if isempty(dir(nameFolder_Figures_perSubj))
                mkdir(nameFolder_Figures_perSubj);
            end

            fprintf('\n%s: %d/%d %s \n', datetime('now'), iFile, numel(nameDir), nameIO)

            % Part 1: templates and metrics
            fprintf('%s: Plotting part 1 (compIV) starts\n', datetime('now'))
            NOMplot_compIV

            pred_test = pred_metrics_allIter{1}; % DV across bins are the same across iterations
            iModelB=iModelB_sim_allCond;
            nParams = numel(namesModelBparams{iModelB});

            % Use the min and max of DV to constrain criterion
            NOMc_lb = min(data_train_allIter{1}.IV);
            NOMc_ub = max(data_train_allIter{1}.IV);
            NOMc_0 = mean([NOMc_lb, NOMc_ub]);
            switch iModelB
                case 1  % FullModel: multi. + additive + shared, criterion
                    params0   = [NOMp1_0,   NOMp2_0,   NOMp3_0, NOMc_0];
                    params_lb = [NOMp1_lb, NOMp2_lb, NOMp3_lb, NOMc_lb];
                    params_ub = [NOMp1_ub, NOMp2_ub, NOMp3_ub, NOMc_ub];

                case 2  % No Nshared: multi. + additive (no shared noise)
                    params0   = [NOMp1_0,   NOMp2_0, NOMc_0];
                    params_lb = [NOMp1_lb, NOMp2_lb, NOMc_lb];
                    params_ub = [NOMp1_ub, NOMp2_ub, NOMc_ub];

                case 3  % No Nmulti: additive + shared (no multi. noise)
                    params0   = [NOMp2_0,   NOMp3_0, NOMc_0];
                    params_lb = [NOMp2_lb, NOMp3_lb, NOMc_lb];
                    params_ub = [NOMp2_ub, NOMp3_ub, NOMc_ub];

                case 4  % No Nadd: multi. + shared (no additive noise)
                    params0   = [NOMp1_0,   NOMp3_0, NOMc_0];
                    params_lb = [NOMp1_lb, NOMp3_lb, NOMc_lb];
                    params_ub = [NOMp1_ub, NOMp3_ub, NOMc_ub];

                otherwise
                    error('OOD_NOM_Trialwise_fitNOM: Unknown iModelB = %d', iModelB);
            end

            % Plot Part 2: estimates
            fprintf('%s: Plotting part 2 (fitNOM) starts\n\n', datetime('now'))
            NOMplot_fitNOM
            figure, text(.5,.5,0, sprintf('True criterion\nIn z-score = %.1f\nIn DV = %.2f\n\n', Cz_true, criterion_DV_true), 'FontSize', 20, 'HorizontalAlignment', 'center')
            saveas(gcf, sprintf('%s/5bTrueCriterion_L%d_A%dB%d.jpg', nameFolder_Figures_perSubj, iLocComb, iModelA, iModelB))
            close all
        end % if flag_plotPerCond

        % Extract the value to compare (truth vs. recovered value)
        switch str_metric
            case 'template' % true vs. estimated template
                valTrue_allIter = truth.template_true(:);
                valData_allIter = data_compIV.template_tmpl_allIter;
                valData_allIter = reshape(valData_allIter, nIter, nORI * nSF);

            case 'metrics' % simulated vs. predicted metrics (pA and pD)
                valTrue_allIter = data_compIV.data_metrics_allIter;
                valData_allIter = data_fitNOM.pred_metrics_allIter; % BUT this is cell!! Need to convert!

            case 'parameters' % true vs. estimated parameters
                valTrue_allIter = truth.(str_params);
                valData_allIter = data_fitNOM.params_est_allIter(:, iParam);
        end

        % Get RMSE and R^2
        % RMSE
        valTrue_allIter_reshape = repmat(valTrue_allIter', nIter, 1); % same shape as valData_allIter: [nORI x nSF] x nIter
        rmse_allIter = sqrt(mean((valTrue_allIter_reshape - valData_allIter).^2, 2));
        [rmse_ave, ~, ~, rmse_sem] = getCI(rmse_allIter, 2, 1);

        % R^2
        R2_allIter = nan;
        if strcmp(str_metric, 'template')
            % Normalize
            valTrue_allIter = valTrue_allIter/max(valTrue_allIter);
            for iIter=1:nIter
                valData_perIter = valData_allIter(iIter, :)/max(valData_allIter(iIter, :));

                % valData_allIter_ = valData_allIter ./ max(valData_allIter, [], 2);
                sst = sumsqr(valData_perIter - mean(valData_perIter));
                sse = sumsqr(valTrue_allIter(:) - valData_perIter(:));
                R2_allIter(iIter) = 1 - sse / sst;
            end
        end

        % Store
        nameIO_allFiles{iFile} = nameIO;
        Cz_allFiles(iFile) = Cz_true;
        lambda_allFiles(iFile) = lambda_true;
        rmse_allFiles(iFile, 1:nIter) = rmse_allIter;
        R2_allFiles(iFile, 1:nIter) = R2_allIter;
        valTrue_allFiles(iFile) = mean(valTrue_allIter); % taking ave is for template
        valData_allFiles(iFile, 1:nIter) = mean(valData_allIter, 2); % taking ave is for template

        % Store estimated Nshared and estimated criterion
        iParam_criterion = 4;
        iParam_Nshared = 3;   % <-- change if Nshared is stored at a different column

        criterion_est_allFiles(iFile, 1:nIter) = data_fitNOM.params_est_allIter(:, iParam_criterion);
        Nshared_est_allFiles(iFile, 1:nIter) = data_fitNOM.params_est_allIter(:, iParam_Nshared);

    end % iFile
end % iModel_sim

rmse_allFiles = rmse_allFiles(:, 1:nIter);
R2_allFiles = R2_allFiles(:, 1:nIter);
valData_allFiles = valData_allFiles(:, 1:nIter);
criterion_est_allFiles = criterion_est_allFiles(:, 1:nIter);
Nshared_est_allFiles = Nshared_est_allFiles(:, 1:nIter);

%% Obtain unique conditions
Cz_unik = unique(Cz_allFiles);
lambda_unik = unique(lambda_allFiles);

%% Make color
absMaxCz = max(abs(Cz_unik));
cmap = nan(numel(Cz_unik), 3);
for iCz = 1:numel(Cz_unik)
    Cz = Cz_unik(iCz);

    if absMaxCz == 0
        strength = 0;
    else
        strength = abs(Cz) / absMaxCz;   % 0 to 1
    end

    lightness = 0.85 - 0.5 * strength;   % near 0 = lighter, far from 0 = darker

    if Cz > 0
        % red scale
        cmap(iCz, :) = [1, lightness, lightness];
    elseif Cz < 0
        % blue scale
        cmap(iCz, :) = [lightness, lightness, 1];
    else
        % exactly zero = gray
        cmap(iCz, :) = [0.5, 0.5, 0.5];
    end
end

%% Plot RMSE/R2 across conditions
figure('Position', [100 100 500 1e3]);

% GoF_mat_all = {rmse_ave_mat, R2_ave_mat};
% GoF_sd_all  = {rmse_sd_mat, R2_sd_mat};
GoF_allFiles = {rmse_allFiles, R2_allFiles};
namesGoF_all    = {'RMSE', 'R^2'};

% Set buffer of x-axis
if numel(lambda_unik) > 1
    dx = mean(diff(lambda_unik));
else
    dx = 0.1;
end

for iMetric = 1:2 % 1=RMSE; 2=R2
    subplot(2,1,iMetric); hold on;

    for iCz = 1:numel(Cz_unik)
        Cz_true = Cz_unik(iCz);
        GoF_allLambda = GoF_allFiles{iMetric}(Cz_allFiles == Cz_true, :);
        [GoF_ave_allLambda , ~, ~, GoF_sem_allLambda] = getCI(GoF_allLambda, 2, 2);
        h = errorbar(lambda_unik, GoF_ave_allLambda, GoF_sem_allLambda, '-o', ...
            'LineWidth', 1.8, ...
            'CapSize', 0, ...
            'MarkerSize', 7, ...
            'Color', cmap(iCz,:), ...
            'MarkerFaceColor', cmap(iCz,:), ...
            'MarkerEdgeColor', 'w', ...
            'DisplayName', sprintf('criterion = %.1f', Cz_unik(iCz)));
    end

    xticks(lambda_unik);
    xlim([min(lambda_unik)-0.5*dx, max(lambda_unik)+0.5*dx]);

    if strcmp(namesGoF_all{iMetric}, 'R^2')
        ylim([.5, 1])
        legend('Location', 'best');
    end
    xlabel('\lambda whitening');
    ylabel(namesGoF_all{iMetric});
    title(namesGoF_all{iMetric});

    set(gca, 'TickDir', 'out');
    box off;
end % iMetric

sgtitle(sprintf('%s recovery error\n nIter=%d', str_metric, nIter))

%% Scatter plot for criterion (true vs. estimated criterion in DV units)

figure('Position', [200 200 600 600]); hold on;

% one point per Cz x lambda condition
valTrue_allFiles   = valTrue_allFiles(:);
[valData_ave, ~, ~, valData_sem_lb, valData_sem_ub] = getCI(valData_allFiles, 1, 2);   % ave across iterations
Cz_allFiles     = Cz_allFiles(:);
lambda_allFiles = lambda_allFiles(:);

% remove invalid rows
% idx_valid = ~isnan(valTrue_allFiles) & ~isnan(valData_ave) & ~isnan(Cz_allFiles) & ~isnan(lambda_allFiles);
% valTrue_allFiles   = valTrue_allFiles(idx_valid);
% valData_ave   = valData_ave(idx_valid);
% Cz_allFiles     = Cz_allFiles(idx_valid);
% lambda_allFiles = lambda_allFiles(idx_valid);

% axis limits
xymin = min([valTrue_allFiles; valData_ave]);
xymax = max([valTrue_allFiles; valData_ave]);
if xymin == xymax
    xymin = xymin - 0.1;
    xymax = xymax + 0.1;
end
dx_text = 0.02 * (xymax - xymin);

% color by Cz sign/magnitude
Cz_unik_plot = unique(Cz_allFiles);
absMaxCz = max(abs(Cz_unik_plot));
cmap_cz = nan(numel(Cz_unik_plot), 3);

for iCz = 1:numel(Cz_unik_plot)
    Cz = Cz_unik_plot(iCz);

    if absMaxCz == 0
        strength = 0;
    else
        strength = abs(Cz) / absMaxCz;
    end

    lightness = 0.85 - 0.5 * strength;

    if Cz > 0
        cmap_cz(iCz,:) = [1, lightness, lightness];
    elseif Cz < 0
        cmap_cz(iCz,:) = [lightness, lightness, 1];
    else
        cmap_cz(iCz,:) = [0.5, 0.5, 0.5];
    end
end

% connect points within each lambda group
lambda_unique = unique(lambda_allFiles);
for iLambda = 1:numel(lambda_unique)
    idxL = lambda_allFiles == lambda_unique(iLambda);

    xL = valTrue_allFiles(idxL);
    yL = valData_ave(idxL);
    ySEM_lb = valData_sem_lb(idxL);
    ySEM_ub = valData_sem_ub(idxL);
    czL = Cz_allFiles(idxL);

    [~, sortIdx] = sort(czL);
    xL = xL(sortIdx);
    yL = yL(sortIdx);

    errorbar(xL, yL, ySEM_lb, ySEM_ub, '.', 'Color', [0.7 0.7 0.7], 'LineWidth', 2.0, 'HandleVisibility', 'off', 'CapSize', 0);
    plot(xL, yL, '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.0, 'HandleVisibility', 'off');

    text(xL(1) - dx_text, yL(1), sprintf('\\lambda = %.1f', lambda_unique(iLambda)), ...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'middle', ...
        'FontSize', 10, ...
        'Color', [0.35 0.35 0.35]);
end

% plot points
for iCz = 1:numel(valTrue_allFiles)
    idxC = find(Cz_unik_plot == Cz_allFiles(iCz), 1);
    scatter(valTrue_allFiles(iCz), valData_ave(iCz), 55, cmap_cz(idxC,:), 'filled', ...
        'MarkerEdgeColor', 'w', 'LineWidth', 0.8, ...
        'HandleVisibility', 'off');
end

% unity line
plot([xymin xymax], [xymin xymax], 'k--', 'LineWidth', 1.5, ...
    'HandleVisibility', 'off');

% legend
hLeg = gobjects(numel(Cz_unik_plot),1);
for iCz = 1:numel(Cz_unik_plot)
    hLeg(iCz) = scatter(nan, nan, 55, cmap_cz(iCz,:), 'filled', ...
        'MarkerEdgeColor', 'w', 'LineWidth', 0.8, ...
        'DisplayName', sprintf('Cz = %.1f', Cz_unik_plot(iCz)));
end

xlabel(sprintf('True %s %s', str_metric, str_params));
ylabel(sprintf('Estimated %s %s', str_metric, str_params));
title(sprintf('Recovery of %s %s\nAve. of %d iterations', str_metric, str_params, nIter));
axis square;
xlim([xymin xymax]);
ylim([xymin xymax]);
box off;
set(gca, 'TickDir', 'out');
legend(hLeg, 'Location', 'best');


%% Trade-off between estimated Nshared and estimated criterion
figure('Position', [200 200 600 600]); hold on;

% one point per file / condition, averaged across iterations
[criterion_est_ave, ~, ~, criterion_est_sem_lb, criterion_est_sem_ub] = getCI(criterion_est_allFiles, 1,2);
[Nshared_est_ave, ~, ~, Nshared_est_sem_lb, Nshared_est_sem_ub] = getCI(Nshared_est_allFiles, 1,2);

% color by criterion sign/magnitude, same rule as before
Cz_unique = unique(Cz_allFiles);
absMaxCz = max(abs(Cz_unique));
cmap_cz = nan(numel(Cz_unique), 3);

for iCz = 1:numel(Cz_unique)
    Cz = Cz_unique(iCz);

    if absMaxCz == 0
        strength = 0;
    else
        strength = abs(Cz) / absMaxCz;
    end

    lightness = 0.85 - 0.5 * strength;

    if Cz > 0
        cmap_cz(iCz, :) = [1, lightness, lightness];
    elseif Cz < 0
        cmap_cz(iCz, :) = [lightness, lightness, 1];
    else
        cmap_cz(iCz, :) = [0.5, 0.5, 0.5];
    end
end

% map each point to its Cz color
cmap_pts = nan(numel(Cz_allFiles), 3);
for iCz = 1:numel(Cz_allFiles)
    idx = find(Cz_unique == Cz_allFiles(iCz), 1);
    cmap_pts(iCz, :) = cmap_cz(idx, :);
end

% set axis limits first so label offsets can use them
xmin = min(criterion_est_ave);
xmax = max(criterion_est_ave);
ymin = min(Nshared_est_ave);
ymax = max(Nshared_est_ave);

if xmin == xmax
    xmin = xmin - 0.1;
    xmax = xmax + 0.1;
end
if ymin == ymax
    ymin = ymin - 0.1;
    ymax = ymax + 0.1;
end

dx_text = 0.02 * (xmax - xmin);

% connect points within each lambda group and label lambda
lambda_unique = unique(lambda_allFiles);
for iLambda = 1:numel(lambda_unique)
    idxL = lambda_allFiles == lambda_unique(iLambda);

    xL = criterion_est_ave(idxL);
    yL = Nshared_est_ave(idxL);
    x_SEMlb = criterion_est_sem_lb(idxL);
    x_SEMub = criterion_est_sem_ub(idxL);
    y_SEMlb = Nshared_est_sem_lb(idxL);
    y_SEMub = Nshared_est_sem_ub(idxL);

    czL = Cz_allFiles(idxL);

    % sort by Cz so lines run from negative to positive criterion
    [~, sortIdx] = sort(czL);
    xL = xL(sortIdx);
    yL = yL(sortIdx);
    x_SEMlb = x_SEMlb(sortIdx);
    x_SEMub = x_SEMub(sortIdx);
    y_SEMlb = y_SEMlb(sortIdx);
    y_SEMub = y_SEMub(sortIdx);

    % errorbar(xL, yL, x_SEMlb, x_SEMub, 'horizontal', '.', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.0, 'HandleVisibility', 'off', 'CapSize', 0);
    % errorbar(xL, yL, y_SEMlb, y_SEMub, 'vertical', '.', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.0, 'HandleVisibility', 'off', 'CapSize', 0);
    plot(xL, yL, '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.0, 'HandleVisibility', 'off');

    % label lambda at the minimum-Cz endpoint
    text(xL(1) - dx_text, yL(1), sprintf('\\lambda = %.1f', lambda_unique(iLambda)), ...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'middle', ...
        'FontSize', 10, ...
        'Color', [0.35 0.35 0.35]);
end

% plot points on top
for iCz = 1:numel(criterion_est_ave)
    scatter(criterion_est_ave(iCz), Nshared_est_ave(iCz), 55, cmap_pts(iCz,:), 'filled', ...
        'MarkerEdgeColor', 'w', 'LineWidth', 0.8, ...
        'HandleVisibility', 'off');
end

% optional reference lines at 0
plot([0 0], [ymin ymax], 'k:', 'LineWidth', 1.0, 'HandleVisibility', 'off');
plot([xmin xmax], [0 0], 'k:', 'LineWidth', 1.0, 'HandleVisibility', 'off');

xlabel('Estimated criterion');
ylabel('Estimated Nshared');
box off;
set(gca, 'TickDir', 'out');

% reset limits after plotting refs
xlim([xmin xmax]);
ylim([ymin ymax]);

% legend for criterion color
hLeg = gobjects(numel(Cz_unique),1);
for iCz = 1:numel(Cz_unique)
    hLeg(iCz) = scatter(nan, nan, 55, cmap_cz(iCz,:), 'filled', ...
        'MarkerEdgeColor', 'w', 'LineWidth', 0.8, ...
        'DisplayName', sprintf('Cz=%.1f', Cz_unique(iCz)));
end
legend(hLeg, 'Location', 'best');
axis square;

% ---------- Partial correlation: Nshared vs criterion, controlling for Cz ----------

% use condition-level means across iterations
criterion_valid = criterion_est_ave;
Nshared_valid   = Nshared_est_ave;
Cz_valid        = Cz_allFiles;
lambda_valid    = lambda_allFiles;

% partial correlation controlling for Cz
[r_partial, p_partial] = partialcorr(criterion_valid, Nshared_valid, Cz_valid);

fprintf('Partial correlation between estimated criterion and estimated Nshared, controlling for Cz:\n');
fprintf('r = %.4f, p = %.4g\n', r_partial, p_partial);

% Cz continuous
T = table(Nshared_valid, criterion_valid, Cz_valid, ...
    'VariableNames', {'Nshared', 'criterion', 'Cz'});
mdl = fitlm(T, 'Nshared ~ criterion + Cz');
disp(mdl);

% Cz categorical
T.Cz_cat = categorical(T.Cz);
mdl = fitlm(T, 'Nshared ~ criterion + Cz_cat');
disp(mdl);

% control for both Cz and lambda
[r_partial2, p_partial2] = partialcorr(criterion_valid, Nshared_valid, [Cz_valid, lambda_valid]);
fprintf('Controlling for Cz and lambda: r = %.2f, p = %.3f\n', r_partial2, p_partial2);

title(sprintf(['Estimated Nshared vs. estimated criterion (Ave. of %d iterations)\n' ...
    'Controlling for Cz and lambda: r = %.2f, p = %.3f\n'], ...
    nIter, r_partial2, p_partial2));