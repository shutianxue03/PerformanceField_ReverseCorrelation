% 
% clc
% close all
% clear all
% warning off
% format compact
% 
% addpath(genpath('Data_OOD'))
% addpath(genpath('fxn_model'))
% addpath(genpath('fxn_simulation'))
% 
% % parpool(12)
% 
% %
% subjList = {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'AS', 'HL', 'FH', 'HA', 'CS','DT'};
% nblocks_allSubj = [200, 240, 210, 240, 220, 220, 215, 220, 205, 210, 205, 165];
% nsubj = length(subjList);
% 
% nB = 1;
% nfiltersOri = 29;
% nfiltersSF = 29;
% nLoc8 = 8;
% nLoc5 = 5;
% itype = 3;
% namesLocComb = {'Fovea', 'LHM', 'UVM', 'RHM', 'LVM', 'HM', 'VM', 'Peri'};
% 
% ntrials_pA= 1e5;
% ni = 1; % number of bootstrapping
% nrep = 2; 
% options = optimoptions('fmincon','MaxIterations', 1e4, 'Display','off');
% 
% %%%%%%%%%%%%%%%%%%%
% %  define and simulate the range
% %%%%%%%%%%%%%%%%%%%
% % Params:
% % 1. internal-external noise ratio of the PRS trials
% % 2. internal-external noise ratio of the ABS trials
% % 3. criterion
% params0 = repmat([.5, .5, 0], nLoc8, 1);
% params_lb = repmat([0, 0, -1], nLoc8, 1);
% params_ub = repmat([1.5, 1.5, 1], nLoc8, 1);
% 
% nparams_model = size(params0, 2);
% 
% %%%%%%%%%%%%%%%%%%%
% % empty comtainers
% %%%%%%%%%%%%%%%%%%%
% dprime_perComb_allSubj = nan(nsubj, nLoc8);
% criterion_perComb_allSubj = dprime_perComb_allSubj;
% pA_PRS_perComb_allSubj = dprime_perComb_allSubj;
% pA_ABS_perComb_allSubj = dprime_perComb_allSubj;
% 
% dprime_est_perLoc_ave_allSubj = nan(nsubj, nLoc8);
% dprime_est_perLoc_std_allSubj = dprime_est_perLoc_ave_allSubj;
% criterion_est_perLoc_ave_allSubj = dprime_est_perLoc_ave_allSubj;
% criterion_est_perLoc_std_allSubj = dprime_est_perLoc_ave_allSubj;
% pA_est_PRS_perLoc_ave_allSubj = dprime_est_perLoc_ave_allSubj;
% pA_est_PRS_perLoc_std_allSubj = dprime_est_perLoc_ave_allSubj;
% pA_est_ABS_perLoc_ave_allSubj = dprime_est_perLoc_ave_allSubj;
% pA_est_ABS_perLoc_std_allSubj = dprime_est_perLoc_ave_allSubj;
% params_est_ave_allSubj = nan(nsubj, nLoc8, nparams_model);
% 
% fprintf('\nNumber of boostraps (nB): %d\nNumber of ORI filters: %d\nNumber of SF filters: %d\n', nB, nfiltersOri, nfiltersSF)
% 
% %
% for isubj = 1:nsubj
%     %%% some files mistakenly contained the wrong info
%     subjList = {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'AS', 'HL', 'FH', 'HA', 'CS','DT'};
%     nblocks_allSubj = [200, 240, 210, 240, 220, 220, 215, 220, 205, 210, 205, 165];
%     nsubj = length(subjList);
%     %%%%%%%%
%     
%     subjName = subjList{isubj};
%     nblocks = nblocks_allSubj(isubj);
%     
%     fprintf('\nSubj name: %s (%d/%d)\nNumber of blocks (nblocks): %d\n', subjName, isubj, nsubj, nblocks)
%     
%     load(sprintf('Data_OOD/%s%d_behavMeas.mat', subjName, nblocks))
%     load(sprintf('Data_OOD/%s%d_raw_%d_%d.mat', subjName, nblocks, nfiltersOri, nfiltersSF))
%     load(sprintf('Data_OOD/%s%d_energy_raw_%d_%d.mat', subjName, nblocks, nfiltersOri, nfiltersSF))
%     
%     nAllTrials = length(energy2D_allT_perComb{1});
%     ntrials = nAllTrials/2;
%     
%     % calculate the internal variable
%     fprintf('Calculating the internal variable... ')
%     IV_PRS_perLoc = cell(nLoc8, 1);
%     IV_ABS_perLoc = cell(nLoc8, 1);
%     
%     figure('Position', [0 0 1800 400])
%     for iLoc = 1:nLoc8
%         % model 1: the implemented template is the 2D kernel
%         template_raw = squeeze(kernels2D_perComb(itype, iLoc, :, :));
%         % model 2: the implemented template is the reconstructed 2D kernel
%         t_min = min(template_raw(:));
%         template_pos = template_raw - t_min; % shift to make all values positive
%         template_recon = mean(template_pos,1) .* mean(template_pos,2);
%         ampFactor = (max(template_raw(:)) - min(template_raw(:)))/( max(template_recon(:)) - min(template_recon(:)));
%         shiftFactor = mean(template_raw(:)) - mean(template_recon(:));
%         template_back = (template_recon + t_min+ shiftFactor)*ampFactor ;
%         
%         %     figure('Position', [0 0 1200 300])
%         %     subplot(1,4,1), imagesc(template_raw), colorbar, axis square, title('raw')
%         %     subplot(1,4,2), imagesc(template_pos), colorbar, axis square, title('shifted up (all pos)')
%         %     subplot(1,4,3), imagesc(template_recon), colorbar, axis square, title('reconstructed')
%         %     subplot(1,4,4), imagesc(template), colorbar, axis square, title('shifted back')
%         
%         template = template_recon; % for YK, using raw generates twp distributions far away from each other
%         % plot the raw/used templates
%         subplot(2,nLoc8, iLoc), imagesc(template_raw), colorbar, axis square, if iLoc == 1, title('2D kernel'), end
%         subplot(2,nLoc8, iLoc+nLoc8), imagesc(template), colorbar, axis square, if iLoc == 1, title('Templated in use'), end
%         
%         % empty containers
%         IV_PRS = nan(ntrials, 1);
%         IV_ABS = IV_PRS ;
%         
%         % get energy profile of all trials
%         e2D = energy2D_allT_perComb{iLoc};
%         
%         % resample trials (feasible? because pA could not be calculated if trials are scrambled)
%         
%         % calculate internal variable
%         parfor itrial = 1:ntrials
%             % method 1: convolve the 2D template and the energy profile
%             e_perT_PRS = squeeze(e2D(itrial, :, :));
%             IV_PRS(itrial) = sum(template(:).*e_perT_PRS(:));
%             e_perT_ABS = squeeze(e2D(itrial+ntrials, :, :));
%             IV_ABS(itrial) = sum(template(:).*e_perT_ABS(:));
%         end % end of itrial
%         IV_PRS_perLoc{iLoc} = IV_PRS;
%         IV_ABS_perLoc{iLoc} = IV_ABS;
%         
%         %         subplot(3,nLoc8, iLoc+nLoc8*2), hold on
%         %         histogram(IV_PRS, 'FaceColor', 'r', 'FaceAlpha', .3)
%         %         histogram(IV_ABS, 'FaceColor', 'b', 'FaceAlpha', .3)
%     end % end of iLoc
%     sgtitle(subjName)
%     
%     fprintf('Done!\n')
%     fileName_perSubj = sprintf('fxn_model/%s.mat', subjName);
%     fileDir_perSubj = dir(fileName_perSubj);
%     
%     if isempty(fileDir_perSubj)
%         %%%%%%%%%%%%%%%%%%%
%         % fit
%         %%%%%%%%%%%%%%%%%%%
%         fprintf('Fitting...')
%         params_est_all = nan([ni, size(params0)]);
%         tic
%         for ii = 1:ni
%             options = optimoptions('fmincon','MaxIterations', 1e4, 'Display','off');
%             fxn_estParams = @(params) fxn_getError(params, dprime_allT_perComb, criterion_allT_perComb, pA_allT_perComb, IV_PRS_perLoc, IV_ABS_perLoc, ntrials_pA);
%             problem_ML = createOptimProblem('fmincon','objective', fxn_estParams,'x0', params0, 'lb', params_lb, 'ub', params_ub, 'options',options);
%             ms_ML = MultiStart('StartPointsToRun', 'bounds', 'Display', 'off', 'UseParallel', 1);
%             
%             params_est = run(ms_ML, problem_ML, nrep);
%             params_est_all(ii, :, :) = params_est;
%             if ~mod(ii, ni/5), fprintf('='), end
%         end
%         fprintf('Done!\n')
%         
%         toc
%         save(fileName_perSubj, 'params_est_all', 'IV_PRS_perLoc', 'IV_ABS_perLoc')
%     else
%         load(fileName_perSubj)
%     end
%     
%     %%%%%%%%%%%%%%%%%%%
%     % make predictions
%     %%%%%%%%%%%%%%%%%%%
%     dprime_est_perLoc_all = nan(ni, nLoc8);
%     criterion_est_perLoc_all = dprime_est_perLoc_all;
%     pA_est_PRS_perLoc_all = dprime_est_perLoc_all;
%     pA_est_ABS_perLoc_all = dprime_est_perLoc_all;
%     
%     for ii = 1:ni
%         params_est = squeeze(params_est_all(ii, :, :));
%         if ii == 1, plotFlag = 1; else, plotFlag = 0; end
%         [dprime_est_perLoc, criterion_est_perLoc, pA_est_PRS_perLoc, pA_est_ABS_perLoc] = fxn_getEstimates(params_est, IV_PRS_perLoc, IV_ABS_perLoc, ntrials_pA, plotFlag);
%         dprime_est_perLoc_all(ii, :) = dprime_est_perLoc;
%         criterion_est_perLoc_all(ii, :) = criterion_est_perLoc;
%         pA_est_PRS_perLoc_all(ii, :) = pA_est_PRS_perLoc;
%         pA_est_ABS_perLoc_all(ii, :) = pA_est_ABS_perLoc;
%     end
%     sgtitle(subjName)
%     
%     %%%%%%%%%%%%%%%%%%%
%     % plot prediction
%     %%%%%%%%%%%%%%%%%%%
%     
%     figure('Position', [0 0 1800 300])
%     % dprime
%     plotPred(1, dprime_allT_perComb, mean(dprime_est_perLoc_all), std(dprime_est_perLoc_all), 'Estimated dprime')
%     legend({'Data' 'Prediction'}, 'Location', 'best')
%     % criterion
%     plotPred(2, criterion_allT_perComb, mean(criterion_est_perLoc_all), std(criterion_est_perLoc_all), 'Estimated criterion', [-1,1])
%     % pA (PRS)
%     plotPred(3, pA_allT_perComb(1:nLoc8, 1), mean(pA_est_PRS_perLoc_all), std(pA_est_PRS_perLoc_all), 'Estimated pA [PRS]', [.5,1])
%     % pA (ABS)
%     plotPred(4, pA_allT_perComb(1:nLoc8, 2), mean(pA_est_ABS_perLoc_all), std(pA_est_ABS_perLoc_all), 'Estimated pA [ABS]', [.5,1])
%     % estimated params
%     params_est_mean = squeeze(nanmean(params_est_all, 1));
%     params_est_std = squeeze(nanstd(params_est_all, [], 1));
%     subplot(1,5,5), hold on
%     errorbar(params_est_mean(:, 1), params_est_std(:, 1),'ro', 'CapSize', 0)
%     errorbar(params_est_mean(:, 2), params_est_std(:, 2),'bo', 'CapSize', 0)
%     yline(1);
%     xticks(1:nLoc8), xticklabels(namesLocComb)
%     xlim([0, nLoc8+1])
%     legend({'PRS', 'ABS'})
%     title('Internal-external noise ratio')
%     sgtitle(subjName)
%     
%     % save results of all subj
%     dprime_perComb_allSubj(isubj, :) = dprime_allT_perComb;
%     criterion_perComb_allSubj(isubj, :) = criterion_allT_perComb;
%     pA_PRS_perComb_allSubj(isubj, :) = pA_allT_perComb(:, 1);
%     pA_ABS_perComb_allSubj(isubj, :) = pA_allT_perComb(:, 2);
%     
%     dprime_est_perLoc_ave_allSubj(isubj, :) = mean(dprime_est_perLoc_all);
%     dprime_est_perLoc_std_allSubj(isubj, :) = std(dprime_est_perLoc_all);
%     criterion_est_perLoc_ave_allSubj(isubj, :) = mean(criterion_est_perLoc_all);
%     criterion_est_perLoc_std_allSubj(isubj, :) = std(criterion_est_perLoc_all);
%     pA_est_PRS_perLoc_ave_allSubj(isubj, :) = mean(pA_est_PRS_perLoc_all);
%     pA_est_PRS_perLoc_std_allSubj(isubj, :) = std(pA_est_PRS_perLoc_all);
%     pA_est_ABS_perLoc_ave_allSubj(isubj, :) = mean(pA_est_ABS_perLoc_all);
%     pA_est_ABS_perLoc_std_allSubj(isubj, :) = std(pA_est_ABS_perLoc_all);
%     params_est_ave_allSubj(isubj, :, :) = params_est_mean;
%     
% end % end of isubj
% 
% %% plot average over all subj
% if nsubj > 1
%     figure('Position', [0 0 1800 300])
%     % dprime
%     plotPred_ave(1, dprime_perComb_allSubj, dprime_est_perLoc_ave_allSubj, 'Estimated dprime')
%     legend({'Data', 'Prediction'}, 'Location', 'best')
%     % criterion
%     plotPred_ave(2, criterion_perComb_allSubj, criterion_est_perLoc_ave_allSubj, 'Estimated criterion', [-1,1])
%     % pA (PRS)
%     plotPred_ave(3, pA_PRS_perComb_allSubj, pA_est_PRS_perLoc_ave_allSubj, 'Estimated pA [PRS]', [.5,1])
%     % pA (ABS)
%     plotPred_ave(4, pA_ABS_perComb_allSubj, pA_est_ABS_perLoc_ave_allSubj, 'Estimated pA [ABS]', [.5,1])
%     % estimated params
%     params_est_allSubj_mean = squeeze(nanmean(params_est_ave_allSubj, 1));
%     params_est_allSubj_std = squeeze(nanstd(params_est_ave_allSubj, [], 1));
%     subplot(1,5,5), hold on
%     errorbar(params_est_allSubj_mean(:, 1), params_est_allSubj_std(:, 1), 'or', 'CapSize', 0)
%     errorbar(params_est_allSubj_mean(:, 2), params_est_allSubj_std(:, 2), 'ob', 'CapSize', 0)
%     yline(1);
%     xticks(1:nLoc8), xticklabels(namesLocComb)
%     xlim([0, nLoc8+1])
%     ylim([.5, 1.5])
%     legend({'PRS', 'ABS'})
%     title('Internal-external noise ratio')
%     
%     sgtitle('All subj')
%     
% end
% 
% % helper function - fxn_getError
% function error = fxn_getError(params, dprime_allT_perComb, criterion_allT_perComb, pA_allT_perComb, IV_PRS_perLoc, IV_ABS_perLoc, ntrials_pA)
% % estimate behav measurements
% [dprime_est_perLoc, criterion_est_perLoc, pA_est_PRS_perLoc, pA_est_ABS_perLoc] = fxn_getEstimates(params, IV_PRS_perLoc, IV_ABS_perLoc, ntrials_pA, 0);
% 
% % calculate errors
% getError = @(a,b) sumsqr((a-b)/(a+b));
% error_dprime = getError(dprime_est_perLoc, dprime_allT_perComb);
% error_criterion = getError(criterion_est_perLoc, criterion_allT_perComb);
% error_pA_PRS = getError(pA_est_PRS_perLoc, pA_allT_perComb(:, 1));
% error_pA_ABS = getError(pA_est_ABS_perLoc, pA_allT_perComb(:, 2));
% error = error_dprime + error_criterion+ error_pA_PRS + error_pA_ABS;
% end
% 
% % helper function - fxn_getEstimates
% function [dprime_est_perLoc, criterion_est_perLoc, pA_est_PRS_perLoc, pA_est_ABS_perLoc] = fxn_getEstimates(params, IV_PRS_perLoc, IV_ABS_perLoc, ntrials_pA, plotFlag)
% SX_normPDF = @(x,mu,sigma) exp(-(x-mu).^2/(2*sigma^2));
% nLoc8 = 8;
% 
% dprime_est_perLoc = nan(nLoc8,1);
% criterion_est_perLoc = dprime_est_perLoc;
% pA_est_PRS_perLoc = dprime_est_perLoc;
% pA_est_ABS_perLoc = dprime_est_perLoc;
% 
% if plotFlag, figure('Position', [0 0 1800 200]), end
% 
% for iLoc = 1:nLoc8
%     % extract params
%     internalN_PRS = params(iLoc, 1); % internal-external noise ratio [PRS]
%     internalN_ABS = params(iLoc, 2);
%     criterion = params(iLoc, 3);
%     %     criterion = 0;
%     
%     IV_PRS_raw = IV_PRS_perLoc{iLoc};
%     IV_ABS_raw = IV_ABS_perLoc{iLoc};
%     
%     % normalize (all data centered at 0 and have SD=1)
%     IV_PRS = (IV_PRS_raw - mean([IV_PRS_raw; IV_ABS_raw]))/std([IV_PRS_raw; IV_ABS_raw]);
%     IV_ABS = (IV_ABS_raw - mean([IV_PRS_raw; IV_ABS_raw]))/std([IV_PRS_raw; IV_ABS_raw]);
%     
%     % calculate the mean/std of the IV distribution
%     mean_IV_PRS = mean(IV_PRS);
%     std_IV_PRS = std(IV_PRS);
%     mean_IV_ABS = mean(IV_ABS);
%     std_IV_ABS = std(IV_ABS);
%     
%     % plot (for checking)
%     if plotFlag, quickplot_sim(iLoc, IV_PRS, IV_ABS, criterion, iLoc), end
%     
%     % get the new std by adding the internal noise
%     std_IV_PRS_n = std_IV_PRS * (1 + internalN_PRS);
%     std_IV_ABS_n = std_IV_ABS * (1 + internalN_ABS);
%     
%     % estimate dprime and criterion
%     pHit = 1-normcdf(criterion, mean_IV_PRS, std_IV_PRS_n);
%     pFA = 1-normcdf(criterion, mean_IV_ABS, std_IV_ABS_n);
%     [dprime_est, criterion_est] = SX_sim06_SDT(pHit, pFA, 1e5); % 1e5 is just an arb number to ensure d'/c not to be Inf
%     dprime_est_perLoc(iLoc) = dprime_est;
%     criterion_est_perLoc(iLoc) = criterion_est;
%     
%     % internal responses (IR) for the first rep [PRS]
%     IR_PRS = randn(2, ntrials_pA) * std_IV_PRS_n + mean_IV_PRS;
%     IR_ABS = randn(2, ntrials_pA) * std_IV_ABS_n + mean_IV_ABS;
%     
%     %     figure, hold on,
%     %     histogram(IR_PRS, 'FaceColor', 'r', 'FaceAlpha', .3, 'EdgeColor', 'none', 'Normalization', 'probability');
%     %     histogram(IR_ABS, 'FaceColor', 'b', 'FaceAlpha', .3, 'EdgeColor', 'none', 'Normalization', 'probability');
%     %     xline(criterion, 'linewidth', 2);
%     %
%     % response (YES or NO)
%     resp_PRS = IR_PRS >= criterion;
%     resp_ABS = IR_ABS <= criterion;
%     % check if the dprime/c from simulation is similar to the theoretical
%     pHit = mean(resp_PRS(:));
%     pFA = 1-mean(resp_ABS(:));
%     dprime_sim =  norminv(pHit) - norminv(pFA);
%     criterion_sim = -(norminv(pHit) + norminv(pFA))/2;
%     %     if (dprime_est - dprime_sim) > eps, error('ALERT: The simulated dprime/c is diff from the theoretical ones.'), end
%     %     if (criterion_est - criterion_sim) > eps, error('ALERT: The simulated dprime/c is diff from the theoretical ones.'), end
%     
%     % estimate pA (by simulation)
%     pA_PRS_sim = mean(resp_PRS(1,:) == resp_PRS(2,:));
%     pA_ABS_sim = mean(resp_ABS(1,:) == resp_ABS(2,:));
%     pA_est_PRS_perLoc(iLoc) = pA_PRS_sim;
%     pA_est_ABS_perLoc(iLoc) = pA_ABS_sim;
%     
%     % estimate pA (by math    
%     %     x = linspace(criterion-4,criterion+4, 1e4);
%     %     phi1 = SX_normPDF(x-(mean_IV_PRS - mean_IV_ABS), 0, sqrt(2)*std_IV_PRS);
%     %     phi2 = normcdf(x, 0, sqrt(2)*internalN_PRS*std_IV_PRS);
%     %     pA_est = mean(phi1 .* (phi2.^2 .* (1-phi2).^2));
%     %     pHit = SX_normPDF(x(x>criterion), mean_iv_prs_perT, std_iv_prs_perT + noiseInt(iLoc, 1));
%     %     pFA = SX_normPDF(x(x<criterion), mean_iv_abs_perT, std_iv_abs_perT + noiseInt(iLoc, 2));
%     %     pCR = SX_normPDF(x(x<criterion), mean_iv_prs_perT, std_iv_prs_perT + noiseInt(iLoc, 1));
%     %     pMiss = SX_normPDF(x(x>criterion), mean_iv_abs_perT, std_iv_abs_perT + noiseInt(iLoc, 2));
%     %     pA_est = sum(pHit.*pFA)+sum(pMiss.*pCR);
%     %     figure, plot(x(x>criterion),pHit, x(x<criterion), pFA, x(x<criterion), pCR, x(x>criterion), pMiss, 'linewidth', 2)
%     %     pA_est_perLoc(iLoc) = pA_est;
%     
%     %     pA_est = (pHit^2 + (1-pFA) ^2 + (1-pHit)^2 + pFA ^2)/4;
%     
%     
% end % end of iLoc
% end
% 
% %% helper fxn - quickplot_sim
% function quickplot_sim(iplot, IV_PRS, IV_ABS, criterion, iLoc)
% %     figure, hold on, histogram(iv_prs, 'FaceColor', 'r', 'FaceAlpha', .3), histogram(iv_abs, 'FaceColor', 'b', 'FaceAlpha', .3)
% namesLocComb = {'Fovea', 'LHM', 'UVM', 'RHM', 'LVM', 'HM', 'VM', 'Peri'};
% 
% subplot(1,8,iplot)
% hold on
% p_prs = histcounts(IV_PRS, 'Normalization', 'probability');
% p_abs = histcounts(IV_ABS, 'Normalization', 'probability');
% histogram(IV_PRS, 'FaceColor', 'r', 'FaceAlpha', .3, 'EdgeColor', 'none', 'Normalization', 'probability');
% histogram(IV_ABS, 'FaceColor', 'b', 'FaceAlpha', .3, 'EdgeColor', 'none', 'Normalization', 'probability');
% xline(mean(IV_PRS), 'r', 'linewidth', 1.5);
% xline(mean(IV_ABS), 'b', 'linewidth', 1.5);
% xline(criterion, 'k', 'linewidth', 2);
% errorbar(mean(IV_PRS), max(p_prs)/2, std(IV_PRS), 'r', 'horizontal', 'CapSize',0, 'linewidth', 1.5)
% errorbar(mean(IV_ABS), max(p_abs)/2, std(IV_ABS), 'b', 'horizontal', 'CapSize',0, 'linewidth', 1.5)
% title(namesLocComb{iLoc})
% end
% 
% %% helper fxn - plotPred
% function plotPred(iplot, data, pred, pred_std, title_, ylim_)
% nLoc8 = 8;
% namesLocComb = {'Fovea', 'LHM', 'UVM', 'RHM', 'LVM', 'HM', 'VM', 'Peri'};
% 
% subplot(1,5, iplot), hold on
% bar(data)
% errorbar(pred, pred_std, 'ko', 'CapSize', 0)
% xticks(1:nLoc8), xticklabels(namesLocComb)
% xlim([0, nLoc8+1])
% if nargin == 6, ylim(ylim_), end
% title(title_)
% end
% 
% %% helper fxn - plotPred_ave
% function plotPred_ave(iplot, data, pred, title_, ylim_)
% nLoc8 = 8;
% namesLocComb = {'Fovea', 'LHM', 'UVM', 'RHM', 'LVM', 'HM', 'VM', 'Peri'};
% 
% data_mean = mean(data);
% data_std= std(data);
% pred_mean = mean(pred);
% pred_std= std(pred);
% 
% subplot(1,5,iplot), hold on
% for iLoc = 1:nLoc8
%     bar(iLoc, data_mean(iLoc), 'FaceColor', 'none')
%     errorbar(iLoc, data_mean(iLoc), data_std(iLoc), 'k.', 'CapSize', 0, 'HandleVisibility', 'off')
%     errorbar(iLoc-.3, pred_mean(iLoc), pred_std(iLoc), 'ko', 'CapSize', 0)
% end
% xticks(1:nLoc8), xticklabels(namesLocComb)
% xlim([0,nLoc8+1])
% if nargin == 5, ylim(ylim_), end
% title(title_)
% end
% 
