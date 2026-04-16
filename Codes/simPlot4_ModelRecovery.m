% %% ============================================================
% % Fig. S8. Model recovery matrix
% % - rows    : fitted models
% % - columns : generating models
% %
% % For each unique simulation condition (excluding generating model B):
% %   1) make heatmap of win rate
% %   2) make heatmap of median delta nLL
% %
% % Then average across conditions (equal weight per condition) and plot:
% %   3) grand-average win-rate heatmap
% %   4) grand-average delta-nLL heatmap
% %
% % Created by Shutian Xue on Apr 4, 2026
% % Revised on Apr 4, 2026
% %% ============================================================
% 
% clear; clc; close all;
% SX_RC1_setting;
% set(0, 'DefaultFigureVisible', 'off')
% 
% %% ---------------- Settings ----------------
% genModelB_all = [1,3,4];   % full, no shared, no multi, no add
% fitModelB_all = [1,3,4];   % same set for fitted models
% iModelA_fit   = 2;     % 1=RC-derived template; 2=true templaye
% tieTol        = 1e-10;
% 
% nGen = numel(genModelB_all);
% nFit = numel(fitModelB_all);
% 
% namesModelGen_all = namesModelB(genModelB_all);
% namesModelFit_all = namesModelB(fitModelB_all);
% 
% % Data folders
% nameFolder_Data = sprintf('%s/Data_Part4', nameFolder_server);
% nameFolder_Data_NOM_Trialwise = sprintf('%s/Data_NOM_Trialwise_%d%d', nameFolder_Data, nORI, nSF);
% 
% % Output folder for figures
% nameFolder_Figures_modelRecovery = sprintf('%s/IO/4ModelRecovery_%d%d', nameFolder_Figures, nORI, nSF);
% if isempty(dir(nameFolder_Figures_modelRecovery))
%     mkdir(nameFolder_Figures_modelRecovery);
% end 
% 
% %% ---------------- Find all generating folders ----------------
% nameDir_all = dir(nameFolder_Data_NOM_Trialwise);
% isSubdir = [nameDir_all.isdir];
% nameDir_all = nameDir_all(isSubdir);
% nameDir_all = {nameDir_all.name};
% nameDir_all = nameDir_all(~ismember(nameDir_all, {'.', '..'}));
% 
% % Keep only folders that look like IO names
% keepIO = startsWith(nameDir_all, 'IO_') & contains(nameDir_all, '_B');
% nameDir_all = nameDir_all(keepIO);
% 
% if isempty(nameDir_all)
%     error('No IO_* folders found in %s', nameFolder_Data_NOM_Trialwise);
% end
% 
% %% ---------------- Parse folder names into condition metadata ----------------
% % folderInfo = struct([]);
% 
% for iFile = 1:numel(nameDir_all)
%     info = parse_nameIO(nameDir_all{iFile});
%     info.folderName = nameDir_all{iFile};
%     folderInfo(iFile) = info;
% end
% 
% % Condition key excludes generating model B
% condKey_all = cell(size(folderInfo));
% for iFile = 1:numel(folderInfo)
%     condKey_all{iFile} = make_condition_key(folderInfo(iFile));
% end
% 
% [condKey_unique, ~, condIdx_all] = unique(condKey_all);
% 
% nCond = numel(condKey_unique);
% fprintf('Found %d unique condition(s) across all folders.\n', nCond);
% 
% %% ---------------- Storage across conditions ----------------
% winRate_cond_all = nan(nFit, nGen, nCond);
% dnLL_med_cond_all = nan(nFit, nGen, nCond);
% dnLL_lb_cond_all  = nan(nFit, nGen, nCond);
% dnLL_ub_cond_all  = nan(nFit, nGen, nCond);
% nSamples_cond_all = nan(nGen, nCond);
% 
% condLabel_all = cell(1, nCond);
% 
% %% ---------------- Loop over conditions ----------------
% for iCond = 1:nCond
% 
%     fprintf('\n=============================================\n');
%     fprintf('Condition %d / %d\n', iCond, nCond);
%     fprintf('=============================================\n');
% 
%     % load nameIO
%     nameFile_cond = matlab.lang.makeValidName(condKey_unique{iCond});
% 
%     % all folders in this condition, across generating models
%     idxCond = find(condIdx_all == iCond);
%     infoCond = folderInfo(idxCond);
% 
%     condLabel_all{iCond} = make_condition_label(infoCond(1));
%     fprintf('%s\n', condLabel_all{iCond});
% 
%     % collect folders by generating model
%     namesFolder_Gen_all = cell(1, nGen);
%     for iGen = 1:nGen
%         thisB = genModelB_all(iGen);
%         keep = arrayfun(@(x) x.B == thisB, infoCond);
%         namesFolder_Gen_all{iGen} = {infoCond(keep).folderName};
% 
%         fprintf('  Generating model: %s (B=%d), found %d folder(s)\n', ...
%             namesModelGen_all{iGen}, thisB, numel(namesFolder_Gen_all{iGen}));
%     end
% 
%     % compute matrices for this condition
%     [winRate_all, dnLL_med_all, dnLL_lb_all, dnLL_ub_all, nSamples] = ...
%         compute_model_recovery_matrix( ...
%             nameFolder_Data_NOM_Trialwise, namesFolder_Gen_all, ...
%             fitModelB_all, iModelA_fit, tieTol);
% 
%     % store
%     winRate_cond_all(:,:,iCond) = winRate_all;
%     dnLL_med_cond_all(:,:,iCond) = dnLL_med_all;
%     dnLL_lb_cond_all(:,:,iCond) = dnLL_lb_all;
%     dnLL_ub_cond_all(:,:,iCond) = dnLL_ub_all;
%     nSamples_cond_all(:,iCond) = nSamples(:);
% 
%     % --------- plot per-condition heatmaps ---------
%     % % 1. Win rate
%     % h1 = figure('Position', [100 100 700 550]); clf;
%     % imagesc(winRate_all);
%     % axis equal tight;
%     % set(gca, ...
%     %     'XTick', 1:nGen, 'XTickLabel', namesModelGen_all, ...
%     %     'YTick', 1:nFit, 'YTickLabel', namesModelFit_all, ...
%     %     'FontSize', 12, 'TickLength', [0 0], ...
%     %     'XAxisLocation', 'top');
%     % title(sprintf('Model recovery: Win rate (%d iterations)\n%s', nSamples(1), condLabel_all{iCond}), 'FontSize', 13, 'FontWeight', 'bold');
%     % xlabel('Generating model');
%     % ylabel('Fitted model');
%     % cb = colorbar; cb.Label.String = 'Win rate';
%     % caxis([0 1]);
%     % hold on;
%     % draw_cell_grid(nGen, nFit);
%     % overlay_numeric_text(winRate_all, '%.2f');
%     % saveas(h1, fullfile(nameFolder_Figures_modelRecovery, sprintf('winRate_%s.png', nameFile_cond)));
%     % close all
%     % 
%     % % 2. Delta nLL
%     % h2 = figure('Position', [850 100 700 550]); clf;
%     % imagesc(dnLL_med_all);
%     % axis equal tight;
%     % set(gca, ...
%     %     'XTick', 1:nGen, 'XTickLabel', namesModelGen_all, ...
%     %     'YTick', 1:nFit, 'YTickLabel', namesModelFit_all, ...
%     %     'FontSize', 12, 'TickLength', [0 0], ...
%     %     'XAxisLocation', 'top');
%     % title(sprintf('Model recovery: median \\Delta nLL (%d iterations)\n%s', nSamples(1), condLabel_all{iCond}), 'FontSize', 13, 'FontWeight', 'bold');
%     % xlabel('Generating model');
%     % ylabel('Fitted model');
%     % cb = colorbar; cb.Label.String = 'Median \Delta nLL';
%     % hold on;
%     % draw_cell_grid(nGen, nFit);
%     % overlay_numeric_text(dnLL_med_all, '%.2f');
%     % saveas(h2, fullfile(nameFolder_Figures_modelRecovery, sprintf('dnLL_%s.png', nameFile_cond)));
%     % close all
% end % iCond
% 
% %% ---------------- Grand average across conditions ----------------
% fprintf('\n\nPlotting averages....\n\n')
% % equal weight per condition: summarize across condition dimension
% [winRate_mean, winRate_lb, winRate_ub] = getCI(winRate_cond_all, 1, 3);
% [dnLL_med_mean, dnLL_mean_lb, dnLL_mean_ub] = getCI(dnLL_med_cond_all, 1, 3);
% 
% % --------- grand-average win rate ---------
% h3 = figure('Position', [100 100 750 600]); clf;
% imagesc(winRate_mean);
% axis equal tight;
% set(gca, ...
%     'XTick', 1:nGen, 'XTickLabel', namesModelGen_all, ...
%     'YTick', 1:nFit, 'YTickLabel', namesModelFit_all, ...
%     'FontSize', 12, 'TickLength', [0 0], ...
%     'XAxisLocation', 'top');
% title(sprintf('Average model recovery across %d conditions: Win rate', nCond), 'FontSize', 14, 'FontWeight', 'bold');
% xlabel('Generating model');
% ylabel('Fitted model');
% cb = colorbar; cb.Label.String = 'Mean win rate across conditions';
% caxis([0 1]);
% hold on;
% draw_cell_grid(nGen, nFit);
% 
% for iFit = 1:nFit
%     for iGen = 1:nGen
%         txt = sprintf('%.2f\n[%.2f, %.2f]', ...
%             winRate_mean(iFit, iGen), ...
%             winRate_lb(iFit, iGen), ...
%             winRate_ub(iFit, iGen));
%         text(iGen, iFit, txt, ...
%             'HorizontalAlignment', 'center', ...
%             'VerticalAlignment', 'middle', ...
%             'FontSize', 10, ...
%             'FontWeight', 'bold', ...
%             'Color', pick_text_color(winRate_mean(iFit, iGen), 0.55));
%     end
% end
% saveas(h3, fullfile(nameFolder_Figures_modelRecovery, 'winRate_aveCond.png'));
% 
% % --------- grand-average delta nLL ---------
% h4 = figure('Position', [900 100 750 600]); clf;
% imagesc(dnLL_med_mean);
% axis equal tight;
% set(gca, ...
%     'XTick', 1:nGen, 'XTickLabel', namesModelGen_all, ...
%     'YTick', 1:nFit, 'YTickLabel', namesModelFit_all, ...
%     'FontSize', 12, 'TickLength', [0 0], ...
%     'XAxisLocation', 'top');
% title(sprintf('Average model recovery across %d conditions: median \\Delta nLL', nCond), 'FontSize', 14, 'FontWeight', 'bold');
% xlabel('Generating model');
% ylabel('Fitted model');
% cb = colorbar; cb.Label.String = 'Mean median \Delta nLL across conditions';
% hold on;
% draw_cell_grid(nGen, nFit);
% 
% for iFit = 1:nFit
%     for iGen = 1:nGen
%         txt = sprintf('%.2f\n[%.2f, %.2f]', ...
%             dnLL_med_mean(iFit, iGen), ...
%             dnLL_mean_lb(iFit, iGen), ...
%             dnLL_mean_ub(iFit, iGen));
%         text(iGen, iFit, txt, ...
%             'HorizontalAlignment', 'center', ...
%             'VerticalAlignment', 'middle', ...
%             'FontSize', 10, ...
%             'FontWeight', 'bold', ...
%             'Color', pick_text_color(dnLL_med_mean(iFit, iGen), median(dnLL_med_mean(:), 'omitnan')));
%     end
% end
% saveas(h4, fullfile(nameFolder_Figures_modelRecovery, 'dnLL_aveCond.png'));
% close all
% 
% %% HELPERS
% function info = parse_nameIO(nameStr)
% % Parse nameIO built in OOD_sim:
% % IO_cN%.0f_cG%.0f_nT%s_Nm%s_Na%s_Ns%s_Cz%.1f_cont%.1f_whiten%.1f_R%d_%d%d_B%d
% 
% expr = ['^IO_cN(?<cN>\d+)_cG(?<cG>\d+)_nT(?<nT>[^_]+)_Nm(?<Nm>[^_]+)_Na(?<Na>[^_]+)_Ns(?<Ns>[^_]+)' ...
%         '_Cz(?<Cz>-?\d+\.?\d*)_cont(?<cont>-?\d+\.?\d*)_whiten(?<whiten>-?\d+\.?\d*)' ...
%         '_R(?<R>\d+)_?(?<nORISF>\d+)_B(?<B>\d+)$'];
% 
% tok = regexp(nameStr, expr, 'names');
% 
% if isempty(tok)
%     error('Could not parse folder name: %s', nameStr);
% end
% 
% info = struct();
% info.cN = tok.cN;              % keep strings for exact matching
% info.cG = tok.cG;
% info.nT = tok.nT;
% info.Nm = tok.Nm;
% info.Na = tok.Na;
% info.Ns = tok.Ns;
% info.Cz = tok.Cz;
% info.cont = tok.cont;
% info.whiten = tok.whiten;
% info.R = tok.R;
% info.nORISF = tok.nORISF;
% info.B = str2double(tok.B);
% end
% 
% function key = make_condition_key(info)
% % Exclude generating model B
% key = sprintf('cN%s_cG%s_nT%s_Nm%s_Na%s_Ns%s_Cz%s_cont%s_whiten%s_R%s_%s', ...
%     info.cN, info.cG, info.nT, info.Nm, info.Na, info.Ns, ...
%     info.Cz, info.cont, info.whiten, info.R, info.nORISF);
% end
% 
% function label = make_condition_label(info)
% label = sprintf(['noise=%s, signal=%s, nTrials=%s, ' ...
%                  'Nmul=%s, Nadd=%s, Nshared=%s, ' ...
%                  'Cz=%s, cont=%s, whiten=%s, R=%s'], ...
%     info.cN, info.cG, info.nT, ...
%     info.Nm, info.Na, info.Ns, ...
%     info.Cz, info.cont, info.whiten, info.R);
% end
% 
% function [winRate_all, dnLL_med_all, dnLL_lb_all, dnLL_ub_all, nSamples] = ...
%     compute_model_recovery_matrix(nameFolder_Data_NOM_Trialwise, namesFolder_Gen_all, ...
%     fitModelB_all, iModelA_fit, tieTol)
% 
% nGen = numel(namesFolder_Gen_all);
% nFit = numel(fitModelB_all);
% 
% winRate_all = nan(nFit, nGen);
% dnLL_med_all = nan(nFit, nGen);
% dnLL_lb_all  = nan(nFit, nGen);
% dnLL_ub_all  = nan(nFit, nGen);
% nSamples = nan(1, nGen);
% 
% for iGen = 1:nGen
%     nameFolder_perGen = namesFolder_Gen_all{iGen};
%     nLL_allSamples = [];
% 
%     for iFolder = 1:numel(nameFolder_perGen)
%         nameFile_current = fullfile(nameFolder_Data_NOM_Trialwise, nameFolder_perGen{iFolder});
%         nLL_thisFolder = cell(1, nFit);
% 
%         for iFit = 1:nFit
%             iModelB_fit = fitModelB_all(iFit);
%             str_file = sprintf('*A%dB%d.mat', iModelA_fit, iModelB_fit);
%             nameDir_perFit = dir(fullfile(nameFile_current, str_file));
% 
%             keep = ~contains({nameDir_perFit.name}, 'min');
%             nameDir_perFit = nameDir_perFit(keep);
% 
%             if isempty(nameDir_perFit)
%                 nLL_thisFolder{iFit} = [];
%                 continue;
%             end
% 
%             data_loaded = load(fullfile(nameDir_perFit(1).folder, nameDir_perFit(1).name), 'nLL_test_allIter');
%             if ~isfield(data_loaded, 'nLL_test_allIter')
%                 nLL_thisFolder{iFit} = [];
%                 continue;
%             end
% 
%             nLL_thisFolder{iFit} = data_loaded.nLL_test_allIter(:);
%         end
% 
%         if any(cellfun(@isempty, nLL_thisFolder))
%             continue;
%         end
% 
%         lenAll = cellfun(@numel, nLL_thisFolder);
%         nUse = min(lenAll);
% 
%         tmp = nan(nUse, nFit);
%         for iFit = 1:nFit
%             tmp(:, iFit) = nLL_thisFolder{iFit}(1:nUse);
%         end
% 
%         nLL_allSamples = [nLL_allSamples; tmp];
%     end
% 
%     badRow = any(~isfinite(nLL_allSamples), 2);
%     nLL_allSamples(badRow, :) = [];
% 
%     nSamples(iGen) = size(nLL_allSamples, 1);
% 
%     if isempty(nLL_allSamples)
%         continue;
%     end
% 
%     nLL_min = min(nLL_allSamples, [], 2);
%     dnLL_allSamples = nLL_allSamples - nLL_min;
% 
%     winnerMat = abs(dnLL_allSamples) <= tieTol;
%     nWinnerPerRow = sum(winnerMat, 2);
%     winnerFrac = winnerMat ./ nWinnerPerRow;
% 
%     for iFit = 1:nFit
%         dnLL_perFit = dnLL_allSamples(:, iFit);
%         winRate_all(iFit, iGen) = mean(winnerFrac(:, iFit), 'omitnan');
%         dnLL_med_all(iFit, iGen) = median(dnLL_perFit, 'omitnan');
%         dnLL_lb_all(iFit, iGen)  = prctile(dnLL_perFit, 2.5);
%         dnLL_ub_all(iFit, iGen)  = prctile(dnLL_perFit, 97.5);
%     end
% end
% end
% 
% function draw_cell_grid(nGen, nFit)
% for iGen = 0.5 : 1 : nGen + 0.5
%     plot([iGen iGen], [0.5 nFit+0.5], 'k-', 'LineWidth', 1);
% end
% for iFit = 0.5 : 1 : nFit + 0.5
%     plot([0.5 nGen+0.5], [iFit iFit], 'k-', 'LineWidth', 1);
% end
% end
% 
% function overlay_numeric_text(M, fmt)
% [nRow, nCol] = size(M);
% for iRow = 1:nRow
%     for iCol = 1:nCol
%         txtColor = pick_text_color(M(iRow, iCol), 0.55);
%         text(iCol, iRow, sprintf(fmt, M(iRow, iCol)), ...
%             'HorizontalAlignment', 'center', ...
%             'VerticalAlignment', 'middle', ...
%             'FontSize', 11, ...
%             'FontWeight', 'bold', ...
%             'Color', txtColor);
%     end
% end
% end
% 
% function c = pick_text_color(val, thresh)
% if isnan(val)
%     c = 'k';
% elseif val > thresh
%     c = 'w';
% else
%     c = 'k';
% end
% end
% 
