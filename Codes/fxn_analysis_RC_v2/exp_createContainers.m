


% NOT IN USE 
 
% switch typeInd
%     case 1 % create allSubj matrix
%         locAll_allSubj = nan(nsubj, nLoc);
%         ntrials_allSubj = nan(nsubj,1);
%         
%         % SDT
%         dprime_allSubj = nan(nsubj, nLoc);
%         criterion_allSubj = dprime_allSubj;
%         pA_allSubj = nan(nsubj, nLoc);
%         RT_log_allSubj = dprime_allSubj;
%         RT_log_mean1_allSubj = nan(nsubj, 2, nLoc);  % 1 = correct trials, 2 = incorrect trials
%         RT_log_mean2_allSubj = nan(nsubj, 2, nLoc);  % 1 = prs trials, 2 = abstrials
%         fft1D_mean_both_allSubj = nan(nsubj, nLoc, 2, stim.aper_psz);
%         
%         % Energy
%         energy2D_allSubj = nan(nsubj, 2, nLoc, nfiltersOri, nfiltersSF); % 2: PRS and ABS
%         rHit_allSubj = nan(nsubj, nLoc, nbins_e);
%         rFA_allSubj = rHit_allSubj;
%         ebin_prs_allSubj = rHit_allSubj;
%         ebin_abs_allSubj = rHit_allSubj;
%         efficiency_allSubj = nan(nsubj, nLoc, 3); % 3 indicates efficincy, d_real and d_ideal
%         
%         % slopes
%         slopes2D_allSubj = nan(nsubj, 2, nLoc, nfiltersOri, nfiltersSF); % 2: PRS and ABS
%         seprblity_slope_allSubj = nan(nsubj, 2, nLoc);
%         
%         % Kernels
%         kernels2D_allSubj = nan(nsubj, ntypes, nLoc, nfiltersOri, nfiltersSF);
%         kernelsSF_allSubj = nan(nsubj, ntypes, nLoc, nfiltersSF);
%         kernelsOri_allSubj = nan(nsubj, ntypes, nLoc, nfiltersOri);
%         seprblity_kernel_allSubj = nan(nsubj, ntypes, nLoc);
%         
%         % Fitting
%         kernelsSF_pred_allSubj = kernelsSF_allSubj;
%         kernelsOri_pred_allSubj = kernelsOri_allSubj;
%         paramsOri_allSubj = nan(nsubj, ntypes, nLoc, length(namesModelParams{1}));
%         paramsSF_allSubj = nan(nsubj, ntypes, nLoc, length(namesModelParams{2}));
%         
%     case 2 % save allSubj data
%         ntrials_allSubj(isubj) = ntrials;
%         
%         % SDT
%         dprime_allSubj(isubj, :) = dprime_perLoc;
%         criterion_allSubj(isubj, :) = criterion_perLoc;
%         pA_allSubj(isubj, :) = mean(pA_perSess_perLoc);
%         RT_log_allSubj(isubj, :) = RT_log_perLoc;
%         RT_log_mean1_allSubj(isubj, :, :) = RT_log_mean1_perLoc; % correct or incorrect trials
%         RT_log_mean2_allSubj(isubj, :, :) = RT_log_mean2_perLoc; % Gabor-prs or abs trials
%         fft1D_mean_both_allSubj(isubj, :, :, :) = fft1D_mean_both;
%         
%         % energy
%         for itype = 1:2
%             if itype==1, itrialStart = 1; itrialEnd = ntrials; else, itrialStart = ntrials+1; itrialEnd = ntrials*2; end
%             for iLoc = 1:nLoc % reorganize
%                 energy2D_allSubj(isubj, itype, iLoc, :, :) = squeeze(mean(energy_perLoc{iLoc}(itrialStart:itrialEnd, :, :),1));
%             end
%         end
%         
%         rHit_allSubj(isubj, :, :) = rHit_tgt;
%         rFA_allSubj(isubj, :, :) = rFA_tgt;
%         ebin_prs_allSubj(isubj, :, :) = ebin_prs_tgt;
%         ebin_abs_allSubj(isubj, :, :) = ebin_abs_tgt;
%         %         efficiency_allSubj(isubj, :, :) = efficiency;
%         
%         % slopes
%         slopes2D_allSubj(isubj, :, :, :, :) = slopes2D_perLoc;
%         seprblity_slope_allSubj(isubj, :, :) = seprblity_slope;
%         
%         % kernels
%         kernels2D_allSubj(isubj, :, :, :, :) = kernels2D;
%         seprblity_kernel_allSubj(isubj, :, :) = seprblity_kernel;
% end