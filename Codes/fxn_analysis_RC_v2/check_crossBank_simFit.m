clc, clear, close all, format compact, warning off
set(0,'DefaultFigureVisible','off')

%% Base settings %%
SX_RC1_setting;

noiseCST=.2; gaborCST=.4; iModelB_sim=1;
% noiseP_true=[0 1e-5 1e-4];

nORI_list = [9, 19, 29]; nBanks = numel(nORI_list);
nCommon = 10; % a common sampling to predict tuning functions from both true and estimated templates
Nmul_true_list = [0, .1, .5];
Nadd_true_list = [0, 5e-5, 1e-4];
Nshared_true_list = [0, 5e-5, 1e-4];
nTrials_list = [4]*1e3;
nC = 3; %check simulate_responses_local() for three levels of z-scored criterion

% For piloting
% nORI_list = [9, 7]; nBanks = numel(nORI_list);
% Nmul_true_list = [0, .1];
% Nadd_true_list = [0, 5e-5];
% Nshared_true_list = [0, 5e-5];
% nTrials_list = [2]*1e3;

%% Loop
% Create placeholders
corr_2D_allComb = cell(numel(Nmul_true_list), numel(Nadd_true_list), numel(Nshared_true_list), nC, numel(nTrials_list));
corr_ORI_allComb = corr_2D_allComb;
corr_SF_allComb = corr_2D_allComb;

for Nmul_true = Nmul_true_list
    for Nadd_true = Nadd_true_list
        for Nshared_true = Nshared_true_list
            for nTrials = nTrials_list

                % Print progress
                time_start = datetime('now');
                fprintf('==============================\n');
                fprintf('\nTime started: %s\nTrue Nmul=%s, Nadd=%s, Nshared=%s, nTrials=%s\n\n', ...
                    time_start, format_num2exp(Nmul_true), format_num2exp(Nadd_true), format_num2exp(Nshared_true), format_num2exp(nTrials))
                fprintf('==============================\n');

                % Nmul_true = noiseP_true(1);
                % Nadd_true = noiseP_true(2);
                % Nshared_true = noiseP_true(3);

                noiseP_true = [Nmul_true, Nadd_true, Nshared_true];

                templateType_true = 1; % 1 = raw;
                IVType_true = 1; % 1 = sum of dot product;
                convolveType_true = 1; % 1 = dot product; 2 = convolution (for fxn_getIV_v3)
                flag_PatchMode_true = 1; % 1 = use target energy ('T')
                flag_permT = 0; %1=permute the input template per trial
                itype_template = 0;
                iFamily_ORI=1; iFamily_SF=2;
                problem_setting = MultiStart('StartPointsToRun', 'bounds','UseParallel', 1, 'Display', 'off');
                nRep=20;
                flag_plot_tuning = 0;

                % Define noise structure
                noise.noiseCST = noiseCST;
                noise.SF_low = 1;
                noise.SF_high = 4;
                noise.SF_low_sampling = noise.SF_low;
                noise.SF_high_sampling= noise.SF_high;
                noise.ppd = 32; % pixels per deg
                noise.psz = 3 * noise.ppd;
                noise.fix_contrast = 1;
                noise.ratio_gaborInTgt= 0.5;
                noise.ratio_base = 0.5;

                % Stim parameters
                ppd = noise.ppd;
                sz_dva = 3;
                sz_pix = sz_dva * ppd;
                signalORI= 90;
                gaborSF = 2;
                gaborSD = 0.8;

                stim.gaborCST = gaborCST;
                stim.psz = sz_dva;
                stim.aper_psz = sz_pix;
                stim.targetOri= signalORI;
                stim.phase = 0;
                stim.gaborSF = gaborSF;
                stim.gaborSD = gaborSD;
                stim.gabor_sz = sz_dva;
                stim.mask = exp_CreateCircularApertureSin(stim);

                limit0to1 = @(x) min(max(x, 0), 1);
                nTrialsPerSess = 100;
                nPairs = nTrials / 2;
                nSess = nTrials / nTrialsPerSess;

                if mod(nTrials, 2) ~= 0
                    error('nTrials must be even (because nPairs = nTrials/2).');
                end

                %% Simulate stimuli once %%
                simStim = generate_raw_trials_local(nTrials, stim, noise);

                %% Main sim x fit loop ==================== %%
                % Set empty placeholders
                xORI_sim = cell(nBanks, 1);
                xSF_sim = xORI_sim;
                xORI_fit = xORI_sim;
                xSF_fit = xORI_sim;

                corr_2D = nan(nBanks, nC, nBanks);
                corr_ORI_Fine = corr_2D;
                corr_SF_Fine = corr_2D;
                margPredORI_true_cell = cell(nBanks, nC, nBanks);
                margPredSF_true_cell = margPredORI_true_cell;
                margPredORI_deriv_cell = margPredORI_true_cell;
                margPredSF_deriv_cell = margPredORI_true_cell;

                for iBank_sim = 1:nBanks
                    nORI_sim = nORI_list(iBank_sim);
                    nSF_sim = nORI_sim;
                    [filtersOri_sim, filtersSF_sim, axis_tuning_sim] = make_bank_from_nORI(nORI_sim, noise);

                    xORI_sim{iBank_sim} = axis_tuning_sim{1};
                    xSF_sim{iBank_sim} = axis_tuning_sim{2};

                    fprintf('\n--------------------------\n');
                    fprintf('Simulating with nORI = %d....\n', nORI_sim);
                    fprintf('----------------------------\n');

                    % --- SIM BANK ---
                    e3D_sim_allT = project_trials_to_energy_local(simStim, filtersOri_sim, filtersSF_sim, fxn_getSigma_SPdomain);

                    template_true_sim = make_true_template_local(stim, filtersOri_sim, filtersSF_sim, fxn_getSigma_SPdomain, templateType_true);
                    template_true_sim = template_true_sim / max(template_true_sim(:)) * 0.2;

                    ORI_bound_sim = [1 nORI_sim];

                    %---------------------%
                    simResp = simulate_responses_local(e3D_sim_allT, template_true_sim, simStim.iPRS_allT, simStim.iPair_allT, noiseP_true, IVType_true, convolveType_true, ORI_bound_sim, flag_permT);
                    %---------------------%

                    fprintf(' — DONE\n')

                    for iC = 1:nC
                        fprintf('\n\n Criterion#%d\n', iC)

                        resp = simResp.resp{iC};
                        criterion_true = simResp.criterion_true{iC};

                        % Define and create folders
                        nameIO = sprintf('IO_cN%.0f_cG%.0f_nT%s_Nm%s_Na%s_Ns%s_c%d_%d%d_B%d', ...
                            noiseCST*100, gaborCST*100, format_num2exp(nTrials), format_num2exp(Nmul_true), format_num2exp(Nadd_true), format_num2exp(Nshared_true), iC, nORI_sim, nSF_sim, iModelB_sim);

                        nameFolder_Figures_perSubj = sprintf('%s/IO/%s', nameFolder_Figures, nameIO);
                        if isempty(dir(nameFolder_Figures_perSubj)), mkdir(nameFolder_Figures_perSubj), end

                        %% Get metrics
                        pYES_sim = mean(resp == 1);
                        pHit_sim = sum((simStim.iPRS_allT == 1) & (resp == 1)) / sum(simStim.iPRS_allT == 1);
                        pFA_sim = sum((simStim.iPRS_allT == 0) & (resp == 1)) / sum(simStim.iPRS_allT == 0);
                        pC_sim = mean( (simStim.iPRS_allT==1 & resp==1) | (simStim.iPRS_allT==0 & resp==0) );

                        %--------------------------------------%
                        [dprime_sim, c_zscore] = SX_sim06_SDT(pHit_sim, pFA_sim);
                        %--------------------------------------%

                        respC_sim = nan(nPairs, 1);
                        for iPairUnik = 1:nPairs
                            respAB = resp(simStim.iPair_allT == iPairUnik);
                            respC_sim(iPairUnik) = (respAB(1) == respAB(2));
                        end
                        pA_sim = mean(respC_sim);

                        %% Plot simulated responses
                        figure('Position', [0 200 600 500]);
                        hold on;
                        histogram(simResp.DV_noisy(simStim.iPRS_allT == 1), 'FaceColor', 'r', 'DisplayName', 'Signal Present', 'normalization', 'probability');
                        histogram(simResp.DV_noisy(simStim.iPRS_allT == 0), 'FaceColor', 'b', 'DisplayName', 'Signal Absent', 'normalization', 'probability');
                        xline(criterion_true, 'LineWidth', 2, 'DisplayName', 'True criterion');
                        xlim([min(simResp.DV_noisy), max(simResp.DV_noisy)]);
                        ylabel('Proportion');
                        legend('show', 'location', 'best');

                        title(sprintf(['CI_{95}=[%.1f, %.1f], Median=%.1f\n[TRUE] gN=%.0f%%, gC=%.0f%%, crit=%.1f, Nmul=%s, Nadd=%s, Nshared=%s' ...
                            '\n[MEASURED] c_z=%.1f, pYES=%.0f%%, pC=%.0f%%, pHit=%.0f%%, pFA=%.0f%%\npA=%.0f%%, corrPass=%.2f'], ...
                            round(quantile(simResp.DV_noisy, [.05, .95, .5]), 1), noiseCST*100, gaborCST*100, criterion_true, format_num2exp(Nmul_true), format_num2exp(Nadd_true), format_num2exp(Nshared_true), ...
                            c_zscore, pYES_sim*100, pC_sim*100, pHit_sim*100, pFA_sim*100, pA_sim*100, simResp.corrPass));

                        saveas(gcf, sprintf('%s/0PerfHist.jpg', nameFolder_Figures_perSubj));
                        close all;

                        if pC_sim <.6,
                            % error('ALERT: Simulated accuracy is %.0f%%, lower than 60%%!', pC_sim*100),
                            flag_lowAcc=1;
                        else
                            flag_lowAcc=0;
                        end
                        fprintf('    Distribution plotted\n')
                        %% Loop through each sampling scale for fitting
                        for iBank_fit = 1:nBanks

                            nORI_fit = nORI_list(iBank_fit);
                            [filtersOri_fit, filtersSF_fit, axis_tuning_fit] = make_bank_from_nORI(nORI_fit, noise);

                            xORI_fit{iBank_fit} = axis_tuning_fit{1};
                            xSF_fit{iBank_fit} = axis_tuning_fit{2};

                            % Print progress report
                            fprintf('\n    Fitting with nORI=%d...', nORI_fit)

                            % --- FIT BANK ---
                            e3D_fit_allT = project_trials_to_energy_local(simStim, filtersOri_fit, filtersSF_fit, fxn_getSigma_SPdomain);

                            template_true_fit = make_true_template_local(stim, filtersOri_fit, filtersSF_fit, fxn_getSigma_SPdomain, templateType_true);
                            template_true_fit = template_true_fit / max(template_true_fit(:)) * 0.2;

                            % derive template using selected trials
                            simStim.iPRS_allT = simStim.raw.iPRS_allT;
                            cst_allT  = repmat(stim.gaborCST, nTrials, 1);

                            switch itype_template
                                case 1
                                    useIdx = (simStim.iPRS_allT == 1);
                                case 2
                                    useIdx = (simStim.iPRS_allT == 0);
                                otherwise
                                    useIdx = true(size(simStim.iPRS_allT));
                            end

                            e3D_fit_sel = e3D_fit_allT(useIdx,:,:);
                            iPRS_sel    = simStim.iPRS_allT(useIdx);
                            resp_sel    = resp(useIdx);
                            cst_sel     = cst_allT(useIdx);

                            if flag_standEnergy
                                e3D_fit_norm = normEnergy(e3D_fit_sel, cst_sel, iPRS_sel);
                            else
                                e3D_fit_norm = e3D_fit_sel;
                            end

                            template_raw_fit = SX_sim07_RC(filtersSF_fit, filtersOri_fit, e3D_fit_norm, resp_sel);
                            template_fit = fxn_getTemplate(template_raw_fit, templateType_true, 0);

                            %% Fit both templates and reconstruct tuning functions on the same sampling grid
                            % True template
                            fitTrue = fit_template_with_tuning_local( ...
                                template_true_sim, axis_tuning_sim, ...
                                iFamily_ORI, iFamily_SF, ...
                                ub_full_all, lb_full_all, ...
                                options_fmin, problem_setting, nRep, flag_plot_tuning, nCommon);

                            % Estimated template
                            fitSim = fit_template_with_tuning_local( ...
                                template_fit, axis_tuning_fit, ...
                                iFamily_ORI, iFamily_SF, ...
                                ub_full_all, lb_full_all, ...
                                options_fmin, problem_setting, nRep, flag_plot_tuning, nCommon);

                            %% Compute and store correlation and pred tuning fxn (for plotting)
                            if flag_lowAcc
                                corr_ORI_Fine(iBank_sim, iC, iBank_fit) = 0;
                                corr_SF_Fine(iBank_sim, iC, iBank_fit) = 0;
                            else
                                if iBank_sim == iBank_fit % calculate the correlation between two templates when their shapes match
                                    corr_2D(iBank_sim, iC, iBank_fit) = corr(template_true_sim(:), template_fit(:));
                                end
                                corr_ORI_Fine(iBank_sim, iC, iBank_fit) = corr(fitTrue.margPred_ORI_fine(:), fitSim.margPred_ORI_fine(:), 'rows', 'complete');
                                corr_SF_Fine(iBank_sim, iC, iBank_fit) = corr(fitTrue.margPred_SF_fine(:), fitSim.margPred_SF_fine(:), 'rows', 'complete');
                            end
                            margPredORI_true_cell{iBank_sim, iC, iBank_fit} = fitTrue.margPred_ORI;
                            margPredSF_true_cell{iBank_sim, iC, iBank_fit}  = fitTrue.margPred_SF;
                            margPredORI_deriv_cell{iBank_sim, iC, iBank_fit} = fitSim.margPred_ORI;
                            margPredSF_deriv_cell{iBank_sim, iC, iBank_fit}  = fitSim.margPred_SF;

                            fprintf('DONE')
                        end % iBank_fit
                    end % iC
                end % iBank_sim

                %% Plot correlation between fitted tuning functions on an identical grid -------------------- %%
                for iC = 1:nC
                    nameCrossBank = sprintf('cN%.0f_cG%.0f_nT%s_Nm%s_Na%s_Ns%s_c%d', ...
                        noiseCST*100, gaborCST*100, format_num2exp(nTrials), format_num2exp(Nmul_true), format_num2exp(Nadd_true), format_num2exp(Nshared_true), iC);

                    %% Plot correlation matrix
                    figure('Position', [200 200 1e3 520]);
                    for iFeature = 1:2
                        switch iFeature
                            case 1
                                corr_Fine = corr_ORI_Fine;
                                str_label = 'nORI';
                                nChannel_list = nORI_list;
                            case 2
                                corr_Fine = corr_SF_Fine;
                                str_label = 'nSF';
                                nChannel_list = nORI_list;
                        end
                        subplot(1,2,iFeature), hold on
                        imagesc(squeeze(corr_Fine(:, iC, :)));
                        axis square; colorbar;
                        xlabel(sprintf('Fit %s', str_label));
                        ylabel(sprintf('Sim %s', str_label));
                        caxis([.5, 1])
                        set(gca, 'XTick', 1:nBanks, 'XTickLabel', string(nChannel_list), 'YTick', 1:nBanks, 'YTickLabel', string(nChannel_list));

                        % Print correlation at the center of each cell
                        for iBank_sim = 1:nBanks
                            for iBank_fit = 1:nBanks
                                % corr between reconstructed tuning fxns
                                val_corr = squeeze(corr_Fine(iBank_sim, iC, iBank_fit));
                                % corr between 2D templates
                                if iBank_sim == iBank_fit
                                    val_corr_2D = squeeze(corr_2D(iBank_sim,iC, iBank_fit));
                                    str_corr = sprintf('%.2f (%.2f)', val_corr, val_corr_2D);
                                else
                                    str_corr = sprintf('%.2f', val_corr);
                                end

                                text(iBank_fit, iBank_sim, sprintf('%s', str_corr), ...
                                    'HorizontalAlignment', 'center', ...
                                    'FontWeight', 'bold', ...
                                    'Color', 'k');
                            end % iBank_fit
                        end % iBank_sim
                    end % iFeature

                    sgtitle(sprintf('Fine-grid reconstructed template correlation\n%s', nameCrossBank));
                    saveas(gcf, sprintf('%s/IO/CrossBank_%s_2D.png', nameFolder_Figures, nameCrossBank))

                    %% Plot tuning functions
                    for iFeature=1:2
                        switch iFeature
                            case 1
                                x_sim = xORI_sim;
                                x_fit = xORI_fit;
                                margPred_true_cell = margPredORI_true_cell;
                                margPred_deriv_cell = margPredORI_deriv_cell;
                                corr_Fine = corr_ORI_Fine;
                                strF = 'ORI';
                            case 2
                                x_sim = xSF_sim;
                                x_fit = xSF_fit;
                                margPred_true_cell = margPredSF_true_cell;
                                margPred_deriv_cell = margPredSF_deriv_cell;
                                corr_Fine = corr_SF_Fine;
                                strF = 'SF';
                        end

                        figure('Position', [200 200 2e3 2e3]);
                        iPlot = 1;
                        for iBank_sim = 1:nBanks
                            for iBank_fit = 1:nBanks
                                subplot(nBanks, nBanks, iPlot), hold on
                                plot(x_sim{iBank_sim}, margPred_true_cell{iBank_sim, iC, iBank_fit}/max(margPred_true_cell{iBank_sim, iC, iBank_fit}), 'o-')
                                plot(x_fit{iBank_fit}, margPred_deriv_cell{iBank_sim, iC, iBank_fit}/max(margPred_deriv_cell{iBank_sim, iC, iBank_fit}), '*-')
                                title(sprintf('nBanks_{sim}=%d, nBanks_{fit}=%d | r=%.2f', ...
                                    numel(x_sim{iBank_sim}), numel(x_fit{iBank_fit}), corr_Fine(iBank_sim, iC, iBank_fit)), ...
                                    'interpreter', 'latex')
                                iPlot = iPlot+1;
                            end
                        end
                        legend({'Simulated', 'Predicted'})
                        saveas(gcf, sprintf('%s/IO/CrossBank_%s_%s.png', nameFolder_Figures, nameCrossBank, strF))
                    end

                    fprintf('\n\n')

                    %% Save r per nTrials, and per parameter combination
                    i_Nmul = find(Nmul_true == Nmul_true_list);
                    i_Nadd = find(Nadd_true == Nadd_true_list);
                    i_Nshared = find(Nshared_true == Nshared_true_list);
                    i_nTrials = find(nTrials == nTrials_list);
                    corr_2D_allComb{i_Nmul, i_Nadd, i_Nshared, iC, i_nTrials} = corr_2D;
                    corr_ORI_allComb{i_Nmul, i_Nadd, i_Nshared, iC, i_nTrials} = corr_ORI_Fine;
                    corr_SF_allComb{i_Nmul, i_Nadd, i_Nshared, iC, i_nTrials} = corr_SF_Fine;
                end % iC

                % Report duration
                time_end = datetime('now')
                elapsed = time_end - time_start;
                fprintf('\n\nDONE (time used: %s)\n\n\n\n', char(elapsed));
            end % nTrials
        end % Nshared_true
    end % Nadd_true
end % Nmul_true

% Save outputs for future analysis
save(sprintf('%s/CrossBank_outputs.mat', nameFolder_Data), 'corr_2D_allComb', 'corr_ORI_allComb', 'corr_SF_allComb')

%% Visualize corr 
% ================= User choice =================
if ~exist('i_nTrials', 'var')
    i_nTrials = 1; % choose which nTrials slice to visualize
end

nBanks   = numel(nORI_list);
nNmul    = numel(Nmul_true_list);
nNadd    = numel(Nadd_true_list);
nNshared = numel(Nshared_true_list);
nC       = size(corr_ORI_allComb, 4); % criterion dimension in cell array

clim_all = [0.5 1.00];   % adjust if needed
flag_show_numbers = 1;    % 1 = print numbers inside mini heatmaps

% ================= Make integrated mini-heatmap figures =================
for iFeature = 1:2

    switch iFeature
        case 1
            corr_allComb = corr_ORI_allComb;
            str_feature = 'ORI';
            str_title_feature = 'Orientation';
        case 2
            corr_allComb = corr_SF_allComb;
            str_feature = 'SF';
            str_title_feature = 'Spatial frequency';
    end

    for i_Nshared = 1:nNshared
        for iC = 1:nC

            figure('Position', [80 80 1200 1000]);
            t = tiledlayout(nBanks, nBanks, 'TileSpacing', 'compact', 'Padding', 'compact');

            for iBank_sim = 1:nBanks
                for iBank_fit = 1:nBanks

                    nexttile; hold on

                    % mini heatmap: rows = Nadd, cols = Nmul
                    miniMat = nan(nNadd, nNmul);

                    for i_Nadd = 1:nNadd
                        for i_Nmul = 1:nNmul

                            % Collect all non-empty values across the redundant cell-dim (4th dim)
                            corr_all = nan(1, size(corr_allComb, 4));

                            for iCcell = 1:size(corr_allComb, 4)
                                tmp_cell = corr_allComb{i_Nmul, i_Nadd, i_Nshared, iCcell, i_nTrials};

                                if isempty(tmp_cell)
                                    continue
                                end

                                % Expected cases:
                                % 1) tmp_cell is nBanks x nC x nBanks
                                % 2) tmp_cell is nBanks x nBanks
                                if ndims(tmp_cell) == 3
                                    corr_all(iCcell) = tmp_cell(iBank_sim, iC, iBank_fit);
                                elseif ismatrix(tmp_cell)
                                    corr_all(iCcell) = tmp_cell(iBank_sim, iBank_fit);
                                end
                            end

                            if all(isnan(corr_all))
                                miniMat(i_Nadd, i_Nmul) = NaN;
                            else
                                [corr_ave, ~, ~] = getCI(corr_all, 2, 2, 0.68, 1); % mean across non-NaN values
                                miniMat(i_Nadd, i_Nmul) = corr_ave;
                            end
                        end
                    end

                    imagesc(miniMat);
                    axis square
                    set(gca, 'YDir', 'normal')
                    caxis(clim_all)

                    set(gca, ...
                        'XTick', 1:nNmul, 'XTickLabel', string(Nmul_true_list), ...
                        'YTick', 1:nNadd, 'YTickLabel', string(Nadd_true_list), ...
                        'FontSize', 8);

                    % only show outer labels to reduce clutter
                    if iBank_sim == nBanks
                        xlabel('Nmul');
                    else
                        set(gca, 'XTickLabel', []);
                    end

                    if iBank_fit == 1
                        ylabel('Nadd');
                    else
                        set(gca, 'YTickLabel', []);
                    end

                    title(sprintf('Sim %d | Fit %d', nORI_list(iBank_sim), nORI_list(iBank_fit)), ...
                        'FontSize', 10);

                    if flag_show_numbers
                        for iRow = 1:nNadd
                            for iCol = 1:nNmul
                                if ~isnan(miniMat(iRow, iCol))
                                    text(iCol, iRow, sprintf('%.2f', miniMat(iRow, iCol)), ...
                                        'HorizontalAlignment', 'center', ...
                                        'FontSize', 7, 'FontWeight', 'bold', 'Color', 'k');
                                end
                            end
                        end
                    end


                end % iBank_fit
            end % iBank_sim

            cb = colorbar;
            cb.Layout.Tile = 'east';
            cb.Label.String = sprintf('Mean correlation (criterion %d)', iC);

            sgtitle(sprintf('%s template recovery | Nshared = %s | criterion %d', ...
                str_title_feature, format_num2exp(Nshared_true_list(i_Nshared)), iC));

            % Save the figure
            nameFolder_IntHeatmap = sprintf('%s/IO/IntegratedHeatmap', nameFolder_Figures);
            if isempty(dir(nameFolder_IntHeatmap)), mkdir(nameFolder_IntHeatmap), end
            saveas(gcf, sprintf('%s/%s_Ns%s_Cz%d.png', ...
                nameFolder_IntHeatmap, str_feature, format_num2exp(Nshared_true_list(i_Nshared)), iC));
        end % iC
    end % i_Nshared
end % iFeature
%% HELPERS ============================

function [filtersOri_all, filtersSF_all, axis_tuning] = make_bank_from_nORI(nORI, noise)
assert(mod(nORI,2)==1, 'nORI must be odd.');

nSF = nORI;

fOri = linspace(0, 90, (nORI+1)/2);
fOri = fOri(2:end);
filtersOri_all = round([-flip(fOri), 0, fOri] + 90);

SF_low = noise.SF_low;
SF_high = 2 / noise.SF_low * 2;
filtersSF_all_log = linspace(log2(SF_low), log2(SF_high), nSF);
filtersSF_all = 2.^filtersSF_all_log;

axis_tuning = cell(1,2);
axis_tuning{1} = filtersOri_all - 90;
axis_tuning{2} = filtersSF_all_log;
end

function raw = generate_raw_trials_local(nTrials, stim, noise)
assert(mod(nTrials,2)==0, 'nTrials must be even.');

nPairs = nTrials/2;
raw.iPRS_allT_OnePass = [ones(nPairs/2,1); zeros(nPairs/2,1)];
raw.iPRS_allT = [raw.iPRS_allT_OnePass; raw.iPRS_allT_OnePass];
iPair_allT = [1:nPairs, 1:nPairs]';

template_gabor_true = exp_CreateGabor(stim, stim.gaborCST);
limit0to1 = @(x) min(max(x, 0), 1);

patch_target_firstPass = cell(nPairs, 1);

parfor iPair = 1:nPairs
    filtered_noise = exp_CreateFilteredNoise(noise);

    if raw.iPRS_allT(iPair) == 1
        gabor = template_gabor_true;
        patch_target = limit0to1(noise.ratio_base + gabor * noise.ratio_gaborInTgt + filtered_noise) .* stim.mask ...
            + noise.ratio_base * (1 - stim.mask);
    else
        patch_target = limit0to1(noise.ratio_base + filtered_noise) .* stim.mask ...
            + noise.ratio_base * (1 - stim.mask);
    end
    patch_target_firstPass{iPair} = patch_target;
end

patch_target_allT = [patch_target_firstPass; patch_target_firstPass];

raw.stim = stim;
raw.noise = noise;
raw.patch_target_allT = patch_target_allT;
raw.raw.iPRS_allT = raw.iPRS_allT;
raw.iPair_allT = iPair_allT;
end

function e3D_target_allT = project_trials_to_energy_local(raw, filtersOri_all, filtersSF_all, fxn_getSigma_SPdomain)
nORI = numel(filtersOri_all);
nSF = numel(filtersSF_all);
[filter_sin, filter_cos] = SX_sim02_setFilters(raw.stim, filtersSF_all, filtersOri_all, fxn_getSigma_SPdomain, 0);

nTrials = numel(raw.patch_target_allT);
e3D_target_allT = nan(nTrials, nORI, nSF);

for it = 1:nTrials
    e3D = SX_RC4_Energy_parfor(raw.stim.mask, raw.patch_target_allT(it), filter_sin, filter_cos);
    e3D_target_allT(it,:,:) = squeeze(e3D);
end

end

function template_true = make_true_template_local(stim, filtersOri_all, filtersSF_all, fxn_getSigma_SPdomain, templateType_true)
[filter_sin, filter_cos] = SX_sim02_setFilters(stim, filtersSF_all, filtersOri_all, fxn_getSigma_SPdomain, 0);
template_gabor_true = exp_CreateGabor(stim, stim.gaborCST);
template_true = SX_RC4_Energy_parfor(stim.mask, {template_gabor_true}, filter_sin, filter_cos);
template_true = squeeze(template_true);
template_true = fxn_getTemplate(template_true, templateType_true, 0);
end

function simResp = simulate_responses_local(e3D_target_allT, template_true, iPRS_allT, iPair_allT, noiseP_true, IVType_true, convolveType_true, ORI_bound, flag_permT)
% IMPORTANT:
% Adjust this wrapper if your updated fxn_getIV_v3 signature differs.
% This assumes your corrected version takes:
%   fxn_getIV_v3(e3D, template, convolveType, IVType, flag_permT, ORI_bound)

DV_clean = fxn_getIV_v3(e3D_target_allT, template_true, convolveType_true, IVType_true, flag_permT, ORI_bound);

% Normalize IV_target
DV_clean = DV_clean/sum(DV_clean(:));

sigmaMul  = noiseP_true(1);
sigmaAadd  = noiseP_true(2);
sigmaShared = noiseP_true(3);

% Private noise SD for each trial/pass:
% combines multiplicative noise (scales with DV) and additive noise
sigmaPriv_allT = sqrt((DV_clean .* sigmaMul).^2 + sigmaAadd^2);

% Draw one independent private-noise sample for every trial/pass
zPriv_allT = randn(size(DV_clean));

% Draw one shared variance for each trial in a pair
% This same value will be used for both passes of the pair
pairIDs = unique(iPair_allT); % Draw one shared z per PAIR, reuse for both passes
zShared_perPair = randn(numel(pairIDs), 1);
% Map each trial to its pair-specific shared fluctuation
[~, idxPair] = ismember(iPair_allT, pairIDs); % Map pair -> z_sh for each trial
zShared_allT = zShared_perPair(idxPair);

% Construct noisy DV
DV_noisy = DV_clean + sigmaShared .* zShared_allT + sigmaPriv_allT .* zPriv_allT;

% calculate correlation between two passes
% rho_t = Nshared_true^2 / (Nshared_true^2 + sigma_private_t^2)
rho = sigmaShared^2/(sigmaShared^2+sigmaPriv_allT.^2);

% Correlation between the noisy DV of two passes (assuming first half and second half of trials are the same trials)
corrPass = corr(DV_noisy(1:numel(DV_noisy)/2), DV_noisy(numel(DV_noisy)/2+1:numel(DV_noisy)));

%%%%
% Titrate a criterion to reach target accuracy pC_titrate
pC_titrate  = .7;
options_fmin = optimoptions('fmincon', 'MaxIterations', 1e4, 'Display', 'off');
%---------------------------%
fxn_loss_pC_ = @(criterion_potential) fxn_loss_pC(criterion_potential, pC_titrate, DV_noisy, iPRS_allT);
%---------------------------%
% criterion_true_mid = fmincon(fxn_loss_pC_, median(DV_noisy), [], [], [], [], min(DV_noisy), max(DV_noisy), [], options_fmin);
% Find the DV-space criterion that yields target c_z (here 0)
c = .5;
c_all = [-c, 0, c];
criterion_true = cell(1, numel(c_all));
resp_allT = cell(1, numel(c_all));
c_z_true = nan(1, numel(c_all));

% Candidate criteria in DV space:
DV_sorted = sort(unique(DV_noisy(:)));
criterion_grid = [-Inf; (DV_sorted(1:end-1) + DV_sorted(2:end))/2; Inf];

for i = 1:numel(c_all)
    target_cz = c_all(i);

    loss_grid = nan(size(criterion_grid));
    c_z_grid = nan(size(criterion_grid));

    for j = 1:numel(criterion_grid)
        c_try = criterion_grid(j);
        [loss_grid(j), c_z_grid(j)] = criterion_loss_cz(c_try, DV_noisy, iPRS_allT, target_cz);
    end

    [~, idx_best] = min(loss_grid);

    criterion_true{i} = criterion_grid(idx_best);
    c_z_true(i) = c_z_grid(idx_best);
    resp_allT{i} = double(DV_noisy > criterion_true{i});
    [~,~, pHit, pFA] = criterion_loss_cz(criterion_grid(idx_best), DV_noisy, iPRS_allT, target_cz);

    % fprintf('Target criterion in z unit = %.2f\n', target_cz);
    % fprintf('Recovered criterion in z unit = %.2f\n', c_z_true(i));
    % fprintf('pHit = %.2f, pFA=%.2f\n\n', pHit, pFA);
end

% figure('Position', [0 200 600 500]);
% hold on;
% histogram(DV_noisy(iPRS_allT == 1), 'FaceColor', 'r', 'DisplayName', 'Signal Present', 'normalization', 'probability');
% histogram(DV_noisy(iPRS_allT == 0), 'FaceColor', 'b', 'DisplayName', 'Signal Absent', 'normalization', 'probability');
% for i = 1:numel(c_all)
%     xline(criterion_true{i}, 'LineWidth', 2);
% end
% xlim([min(DV_noisy), max(DV_noisy)]);
% ylabel('Proportion');
% legend('show', 'location', 'best');

simResp.IV_target = DV_clean;
simResp.corrPass = corrPass;
simResp.DV_noisy = DV_noisy;
simResp.resp = resp_allT;
simResp.criterion_true = criterion_true;
end


function fitOut = fit_template_with_tuning_local(template_InUse, axis_tuning, iFamily_ORI, iFamily_SF, ub_full_all, lb_full_all, options_fmin, problem_setting, nRep, flag_plot_tuning, nFine)

margORI = mean(template_InUse, 2)';
margSF  = mean(template_InUse, 1);

% ----- ORI fit -----
xORI = axis_tuning{1};
fxn_tuningLoss_ORI = @(param_est) sum((margORI - predSFkernel(xORI, iFamily_ORI, param_est, flag_plot_tuning)).^2);

problem_ORI = createOptimProblem('fmincon', ...
    'objective', fxn_tuningLoss_ORI, ...
    'x0', (ub_full_all{iFamily_ORI}+lb_full_all{iFamily_ORI})/2, ...
    'lb', lb_full_all{iFamily_ORI}, ...
    'ub', ub_full_all{iFamily_ORI}, ...
    'options', options_fmin);

margParams_ORI = run(problem_setting, problem_ORI, nRep);
margPred_ORI = predSFkernel(xORI, iFamily_ORI, margParams_ORI, flag_plot_tuning);
margR2_ORI = 1 - sumsqr(margORI - margPred_ORI) / sumsqr(margORI - mean(margORI));

% ----- SF fit -----
xSF_log = axis_tuning{2};
xSF_lin = 2.^xSF_log;
fxn_tuningLoss_SF = @(param_est) sum((margSF - predSFkernel(xSF_lin, iFamily_SF, param_est, flag_plot_tuning)).^2);

problem_SF = createOptimProblem('fmincon', ...
    'objective', fxn_tuningLoss_SF, ...
    'x0', (ub_full_all{iFamily_SF}+lb_full_all{iFamily_SF})/2, ...
    'lb', lb_full_all{iFamily_SF}, ...
    'ub', ub_full_all{iFamily_SF}, ...
    'options', options_fmin);

margParams_SF = run(problem_setting, problem_SF, nRep);
margPred_SF = predSFkernel(xSF_lin, iFamily_SF, margParams_SF, flag_plot_tuning);
margR2_SF = 1 - sumsqr(margSF - margPred_SF) / sumsqr(margSF - mean(margSF));

% ----- Fine reconstruction -----
xORI_fine = linspace(min(xORI), max(xORI), nFine);
xSF_log_fine = linspace(min(xSF_log), max(xSF_log), nFine);
xSF_lin_fine = 2.^xSF_log_fine;

margPred_ORI_fine = predSFkernel(xORI_fine, iFamily_ORI, margParams_ORI, flag_plot_tuning);
margPred_SF_fine  = predSFkernel(xSF_lin_fine, iFamily_SF, margParams_SF, flag_plot_tuning);

templateFine = margPred_ORI_fine(:) * margPred_SF_fine(:)';

fitOut.xORI_fine = xORI_fine;
fitOut.xSF_log_fine = xSF_log_fine;
fitOut.margParams_ORI = margParams_ORI;
fitOut.margParams_SF  = margParams_SF;
fitOut.margPred_ORI = margPred_ORI;
fitOut.margPred_SF = margPred_SF;
fitOut.margPred_ORI_fine = margPred_ORI_fine;
fitOut.margPred_SF_fine = margPred_SF_fine;
fitOut.margR2_ORI = margR2_ORI;
fitOut.margR2_SF  = margR2_SF;
fitOut.templateFine = templateFine;
end

function [loss, c_z, pHit, pFA] = criterion_loss_cz(k, DV_noisy, iPRS_allT, target_cz)

nPRS = sum(iPRS_allT == 1);
nABS = sum(iPRS_allT == 0);
eps_ = 0.5 / min(nPRS, nABS);

resp = double(DV_noisy > k);

pHit = sum(resp == 1 & iPRS_allT == 1) / nPRS;
pFA  = sum(resp == 1 & iPRS_allT == 0) / nABS;

% clip rates to avoid infinities
pHit = min(max(pHit, eps_), 1 - eps_);
pFA  = min(max(pFA,  eps_), 1 - eps_);

c_z = -0.5 * (norminv(pHit) + norminv(pFA));

loss = (c_z - target_cz).^2;
end