% 
% params_est_allSubj = nan(nsubj, nLoc, nparams_model);
% 
% %% plot IDVD data
% for isubj = 1:nsubj
%     
%     subjName = subjList{isubj};
%     nblocks = nblocks_allSubj(isubj);
%     fprintf('\nSubj name: %s (%d/%d)\nNumber of blocks (nblocks): %d\n', subjName, isubj, nsubj, nblocks)
%     
%     % load data
%     if renameMode
%         fileName_perSubj = sprintf('Data_OOD/%s_model.mat', subjName);
%     else
%         fileName_perSubj = sprintf('Data_OOD/%s_model%d.mat', subjName, modelVersion);
%     end
%     load(fileName_perSubj)
%     
%     %% plot the template
%     figure('Position', [1000 200 225*nLoc 300])
%     for iLoc = 1:nLoc
%         subplot(2, nLoc, iLoc), RCplot_2Dkernel(template_raw{iLoc}.', []), if iLoc == 1, title('raw'), end
%         subplot(2, nLoc, iLoc+nLoc), RCplot_2Dkernel(template_recon{iLoc}.', []), colorbar, axis square, if iLoc == 1,title('reconstructed'), end
% %         subplot(3, nLoc, iLoc+nLoc*2), imagesc(template_back{iLoc}), colorbar, axis square, if iLoc == 1,title('shifted back'), end
%     end
%     sgtitle(subjName)
%     
%     %% plot histograms
%     figure('Position', [1000 200 225*nLoc 200])
%     for iLoc = 1:nLoc
%         IV_PRS = data.IV_PRS{iLoc};
%         IV_ABS = data.IV_ABS{iLoc};
%         IV_PRS_norm = (IV_PRS - mean([IV_PRS; IV_ABS]))/std([IV_PRS; IV_ABS]);
%         IV_ABS_norm = (IV_ABS - mean([IV_PRS; IV_ABS]))/std([IV_PRS; IV_ABS]);
%         
%         ModelPlot_hist(iLoc, nLoc, IV_PRS_norm, IV_ABS_norm, params_est_all)
%     end
%     sgtitle(sprintf('Histogram of internal variable - %s', subjName))
%     
% 
%     %% plot prediction per IDVD
%     fxn_plotPred(iLoc_all, 1, data, pred_all, params_est_all, subjName) % do not change nsubj=1!!
%     
%     %%%%%%%%%%%%%%%%%%%
%     % save results of all subj
%     %%%%%%%%%%%%%%%%%%%
%     
%     data_allSubj.dprime(isubj, :) = data.dprime;
%     data_allSubj.criterion(isubj, :) = data.criterion;
%     data_allSubj.pA(isubj, :) = data.pA;
%     data_allSubj.pA_PRS(isubj, :) = data.pA_PRS;
%     data_allSubj.pA_ABS(isubj, :) = data.pA_ABS;
%     data_allSubj.pC(isubj, :) = data.pC;
%     data_allSubj.pHit(isubj, :) = data.pHit;
%     data_allSubj.pFA(isubj, :) = data.pFA;
%     
%     getAVE = @(data) squeeze(median(data, 1));
%     pred_med_allSubj.dprime(isubj, :) = getAVE(pred_all.dprime);
%     pred_med_allSubj.criterion(isubj, :) = getAVE(pred_all.criterion);
%     pred_med_allSubj.pA(isubj, :) = getAVE(pred_all.pA);
%     pred_med_allSubj.pA_PRS(isubj, :) = getAVE(pred_all.pA_PRS);
%     pred_med_allSubj.pA_ABS(isubj, :) = getAVE(pred_all.pA_ABS);
%     pred_med_allSubj.pC(isubj, :) = getAVE(pred_all.pC);
%     pred_med_allSubj.pHit(isubj, :) = getAVE(pred_all.pHit);
%     pred_med_allSubj.pFA(isubj, :) = getAVE(pred_all.pFA);
%     
%     params_est_allSubj(isubj, :, :) = getAVE(params_est_all); % nsubj x ni x nLoc x nparams_model
%     
%     %%%%%%%%%%%%%%%%%%%
%     % rename file
%     %%%%%%%%%%%%%%%%%%%
%     if renameMode
%         fileName_perSubj_new = sprintf('Data_OOD/%s_model%d.mat', subjName, modelVersion);
%         movefile(fileName_perSubj, fileName_perSubj_new);
%     end
%     
% end % end of isubj
% 
% 
% %% plot average over all subj
% if nsubj > 1
%     fileName_allSubj = sprintf('Data_OOD/n%d_model.mat', nsubj);
%     save(fileName_allSubj, 'nsubj', 'nLoc', 'data_allSubj', 'pred_med_allSubj', 'params_est_allSubj')
%     fxn_plotPred(iLoc_all, nsubj, data_allSubj, pred_med_allSubj, params_est_allSubj, sprintf('[Model %d] n = %d', modelVersion, nsubj))
%     
%     %% plot the correlation between pA and IE
%     
%     figure
%     for iiLoc = 1:nLoc
%         iLoc = iLoc_all(iiLoc);
%         for itype = 1:2
%             subplot(nLoc, 2, 2*(iiLoc-1)+itype), hold on
%             pA = squeeze(pA_allT_allSubj(:, iiLoc, itype));
%             alpha = params_est_allSubj(:, iiLoc, itype);
%             % correlation
%             [r, p] = corr(pA, alpha);
%             
%             % linear regression
%             [beta, stats] = polyfit(pA, alpha,1);
%             slope = beta(1);
%             yfit = polyval(beta, pA);
%             plot(pA, yfit, '-', 'color', colors_comb(iLoc, :), 'handlevisibility', 'off')
%             % plot
%             plot(pA, alpha, 'o', 'color', colors_comb(iLoc, :))
%             xlabel(sprintf('pA %s', namesType{itype}))
%             ylabel(sprintf('alpha %s', namesType{itype}))
%             xlim([.5, .9])
%             if itype == 1, ylim([0,2]), end
%             title(sprintf('[%s] %s, r=%.2f, p=%.2f', namesType{itype}, namesLocComb{iLoc}, r, p))
%         end
%     end
%     
% end
