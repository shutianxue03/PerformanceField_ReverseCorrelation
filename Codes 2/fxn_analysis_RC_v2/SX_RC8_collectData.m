
nlines = 2;
switch dataMode
    case 1
        %% create containers for GROUP data [_allSubj]
        %         nAllTrials_allSubj = nan(nsubj, 1);
        nm10 = 10; % (1) dprime (2) criterion (3) RT (4) CST (5-7) pA _all/PRS/ABS (8-10) pC/pHit/pFA
        data_aveBlk_allSubj = nan(nsubj, nm10, nLoc8);
        data_semBlk_allSubj = data_aveBlk_allSubj;
        dprime_allT_allSubj = nan(nsubj, nLoc8);
        criterion_allT_allSubj = dprime_allT_allSubj;
        pA_allT_allSubj = nan(nsubj, nLoc8, ntypes);
        pC_allT_allSubj = pA_allT_allSubj;
        
        % slopes
        pC_tgt_allSubj =  nan([nsubj, size(pC_tgt_allB)]);
        ebin_tgt_allSubj = pC_tgt_allSubj;
        pC_tgt_norm_allSubj = pC_tgt_allSubj;
        ebin_tgt_norm_allSubj = pC_tgt_allSubj;
        
        % kernels
        kernels2D_allSubj = nan([nsubj, size(kernels2D_allB)]);
        intercept2D_allSubj = nan([nsubj, size(intercept2D_allB)]);
        margORI_allSubj = nan([nsubj, size(margORI_allB)]);
        margSF_allSubj = nan([nsubj, size(margSF_allB)]);
        R2_2D_allSubj = nan([nsubj, size(R2_2D_allB)]);
        R2_Tjur_allSubj = R2_2D_allSubj;
        pValues_allSubj = nan([nsubj, size(pValues_allB)]);
        pCat_allSubj = R2_2D_allSubj;
        sep_allSubj =  nan([nsubj, size(sep_allB)]);
        
        % fitting
        % ORI - kernel fitting
        margPred_ORI_allSubj = nan([nsubj, size(margPred_ORI_allB)]);
        margParams_ORI_allSubj = nan([nsubj, size(margParams_ORI_allB)]);
        margR2_ORI_allSubj = nan([nsubj, size(margR2_ORI_allB)]);
        % SF - kernel fitting
        margPred_SF_allSubj = nan([nsubj, size(margPred_SF_allB)]);
        margParams_SF_allSubj = nan([nsubj, size(margParams_SF_allB)]);
        margR2_SF_allSubj = nan([nsubj, size(margR2_SF_allB)]);
        
    case 2
        %% create containers for IDVD data [_allB]
        % % linear regression slopes - behavioral
        pC_tgt_allB = nan(nB, ntypes, nLoc8, nbins_e);
        ebin_tgt_allB = pC_tgt_allB;
        pC_tgt_allB_norm = pC_tgt_allB;
        ebin_tgt_allB_norm = pC_tgt_allB;
        
        % kernel
        kernels2D_allB = nan(nB, ntypes, nLoc8, nfiltersOri, nfiltersSF); % ntypes, nLoc8, nfiltersOri, nfiltersSF
        intercept2D_allB = kernels2D_allB;
        margORI_allB = nan(nB, ntypes, nfiltersOri, nLoc8);
        margSF_allB = nan(nB, ntypes, nfiltersSF, nLoc8);
        R2_2D_allB = nan(nB, ntypes, nLoc8, nfiltersOri, nfiltersSF);
        R2_Tjur_allB = R2_2D_allB;
        pValues_allB = nan(nB, ntypes, nLoc8, nfiltersOri, nfiltersSF, 2); % ntypes, nLoc8, nfiltersOri, nfiltersSF, 2; % 2 = p of intercept+p of slope
        pCat_allB = R2_2D_allB;
        sep_allB = nan(nB, ntypes, nLoc8); % ntypes, nLoc8
        
        % ORI-kernel fitting
        margPred_ORI_allB = nan(nB, ntypes, ncomb6, nlines, nfiltersOri*nitp);
        margParams_ORI_allB = nan(nB, ntypes, ncomb6, nlines, nparamsOri);
        margR2_ORI_allB = nan(nB, ntypes, ncomb6, nlines);
        
        % SF - kernel fitting
        margPred_SF_allB = nan(nB, ntypes, ncomb6, nlines, nfiltersSF*nitp);
        margParams_SF_allB = nan(nB, ntypes, ncomb6, nlines, nparamsSF);
        margR2_SF_allB = margR2_ORI_allB;
        
    case 3
        %% save IDVD data [_allB]
        margOri_perComb = cat(3, margOri_perComb{:});
        margSF_perComb = cat(3, margSF_perComb{:});
        
        % Slope (3D, ntypes x nLoc8 x nbins )
        pC_tgt_allB(iB, :, :, :) = pC_tgt_perComb;
        ebin_tgt_allB(iB, :, :, :) = ebin_tgt_perComb;
        pC_tgt_norm_allB(iB, :, :, :) = pC_tgt_norm_perComb;
        ebin_tgt_norm_allB(iB, :, :, :) = ebin_tgt_norm_perComb;
        % Kernel (4D, ntypes x nLoc8 x nfiltersOri x nfiltersSF, except pValues (5D, last dim is 2) and sep (2D, ntypes x nLoc8))
        kernels2D_allB(iB, :, :, :, :) = kernels2D_perComb;
        intercept2D_allB(iB, :, :, :, :) = intercept2D_perComb;
        R2_2D_allB(iB, :, :, :, :) = R2_2D_perComb;
        R2_Tjur_allB(iB, :, :, :, :) = R2_Tjur_perComb;
        pValues_allB(iB, :, :, :, :, :) = pValues_perComb; % 2 = p of intercept+p of slope
        pCat_allB(iB, :, :, :, :) = pCat_perComb;
        sep_allB(iB, :, :) = sep_perComb; % ntypes x nLoc8
        % MARG (3D, ntypes x nfilters x nLoc8)
        margORI_allB(iB, :, :, :) = margOri_perComb;
        margSF_allB(iB, :, :, :) = margSF_perComb;
        % kernel fitting (ORI)
        margPred_ORI_allB(iB, :, :, :, :) = margPred_ORI_perComb;        % 4D, ntypes x ncomb6 x nlines x nfilters_itp
        margParams_ORI_allB(iB, :, :, :, :) = margParams_ORI_perComb; % 4D, ntypes x ncomb6 x nlines x nparams
        margR2_ORI_allB(iB, :, :, :) = margR2_ORI_perComb;                 % 3D, ntypes x ncomb6 x nlines
        % kernel fitting (SF)
        margPred_SF_allB(iB, :, :, :, :) = margPred_SF_perComb;
        margParams_SF_allB(iB, :, :, :, :) = margParams_SF_perComb;
        margR2_SF_allB(iB, :, :, :) = margR2_SF_perComb;
        
    case 4
        %% save all subj
        % reorder the location (to select the dominant eye side for HM)
        if ~isnan(eyeDom)
            % slope
            pC_tgt_allB = reorderLoc(pC_tgt_allB, 1, eyeDom(isubj));
            ebin_tgt_allB = reorderLoc(ebin_tgt_allB, 1, eyeDom(isubj));
            pC_tgt_norm_allB = reorderLoc(pC_tgt_norm_allB, 1, eyeDom(isubj));
            ebin_tgt_norm_allB = reorderLoc(ebin_tgt_norm_allB, 1, eyeDom(isubj));
            % kernel
            kernels2D_allB = reorderLoc(kernels2D_allB, 2, eyeDom(isubj));
            intercept2D_allB = reorderLoc(intercept2D_allB, 2, eyeDom(isubj));
            R2_2D_allB = reorderLoc(R2_2D_allB, 2, eyeDom(isubj));
            R2_Tjur_allB = reorderLoc(R2_Tjur_allB, 2, eyeDom(isubj));
            pValues_allB = reorderLoc(pValues_allB, 3, eyeDom(isubj));
            pCat_allB = reorderLoc(pCat_allB, 2, eyeDom(isubj));
            sep_allB = reorderLoc(sep_allB, 4, eyeDom(isubj));
            % MARG
            margORI_allB = reorderLoc(margORI_allB, 5, eyeDom(isubj));
            margSF_allB = reorderLoc(margSF_allB, 5, eyeDom(isubj));
            % kernel fitting (ORI)
            margPred_ORI_allB = reorderLoc(margPred_ORI_allB, 6, eyeDom(isubj));
            margParams_ORI_allB = reorderLoc(margParams_ORI_allB, 6, eyeDom(isubj));
            margR2_ORI_allB = reorderLoc(margR2_ORI_allB, 7, eyeDom(isubj));
            % kernel fitting (SF)
            margPred_SF_allB = reorderLoc(margPred_SF_allB, 6, eyeDom(isubj));
            margParams_SF_allB = reorderLoc(margParams_SF_allB, 6, eyeDom(isubj));
            margR2_SF_allB = reorderLoc(margR2_SF_allB, 7, eyeDom(isubj));
        end
        
        % behav
        data_aveBlk_allSubj(isubj, :, :) = data_aveBlk_perComb;
        data_semBlk_allSubj(isubj, :, :) = data_semBlk_perComb;
        dprime_allT_allSubj(isubj, :) = dprime_allT_perComb;
        criterion_allT_allSubj(isubj, :) = criterion_allT_perComb;
        pA_allT_allSubj(isubj, :, :) = pA_allT_perComb;
        pC_allT_allSubj(isubj, :, :) = pC_allT_perComb;
        
        % slope
        pC_tgt_allSubj(isubj, :, :, :, :) = pC_tgt_allB;
        ebin_tgt_allSubj(isubj, :, :, :, :) = ebin_tgt_allB;
        pC_tgt_norm_allSubj(isubj, :, :, :, :) = pC_tgt_norm_allB;
        ebin_tgt_norm_allSubj(isubj, :, :, :, :) = ebin_tgt_norm_allB;
        % kernel
        kernels2D_allSubj(isubj, :, :, :, :, :) = kernels2D_allB;
        intercept2D_allSubj(isubj, :, :, :, :, :) = intercept2D_allB;
        R2_2D_allSubj(isubj, :, :, :, :, :) = R2_2D_allB;
        R2_Tjur_allSubj(isubj, :, :, :, :, :) = R2_Tjur_allB;
        pValues_allSubj(isubj, :, :, :, :, :, :) = pValues_allB;
        pCat_allSubj(isubj, :, :, :, :, :) = pCat_allB;
        sep_allSubj(isubj, :, :, :) = sep_allB;
        % MARG
        margORI_allSubj(isubj, :, :, :, :) = margORI_allB;
        margSF_allSubj(isubj, :, :, :, :) = margSF_allB;
        % ORI - kernel fitting
        margPred_ORI_allSubj(isubj, :, :, :, :, :) = margPred_ORI_allB;
        margParams_ORI_allSubj(isubj, :, :, :, :, :) = margParams_ORI_allB;
        margR2_ORI_allSubj(isubj, :, :, :, :) = margR2_ORI_allB;
        % SF - kernel fitting
        margPred_SF_allSubj(isubj, :, :, :, :, :) = margPred_SF_allB;
        margParams_SF_allSubj(isubj, :, :, :, :, :) = margParams_SF_allB;
        margR2_SF_allSubj(isubj, :, :, :, :) = margR2_SF_allB;
    case 5
        %% compile boostrapped otput of each oatch
        % Slope (3D, ntypes x nLoc8 x nbins )
        pC_tgt_allBB(iB_start:iB_end, :, :, :) = pC_tgt_allB;
        ebin_tgt_allBB(iB_start:iB_end, :, :, :) = ebin_tgt_allB;
        pC_tgt_norm_allBB(iB_start:iB_end, :, :, :) = pC_tgt_norm_allB;
        ebin_tgt_norm_allBB(iB_start:iB_end, :, :, :) = ebin_tgt_norm_allB;
        % Kernel (4D, ntypes x nLoc8 x nfiltersOri x nfiltersSF, except pValues (5D, last dim is 2) and sep (2D, ntypes x nLoc8))
        kernels2D_allBB(iB_start:iB_end, :, :, :, :) = kernels2D_allB;
        intercept2D_allBB(iB_start:iB_end, :, :, :, :) = intercept2D_allB;
        R2_2D_allBB(iB_start:iB_end, :, :, :, :) = R2_2D_allB;
        R2_Tjur_allBB(iB_start:iB_end, :, :, :, :) = R2_Tjur_allB;
        pValues_allBB(iB_start:iB_end, :, :, :, :, :) = pValues_allB; % 2 = p of intercept+p of slope
        pCat_allBB(iB_start:iB_end, :, :, :, :) = pCat_allB;
        sep_allBB(iB_start:iB_end, :, :) = sep_allB; % ntypes x nLoc8
        % MARG (3D, ntypes x nfilters x nLoc8)
        margORI_allBB(iB_start:iB_end, :, :, :) = margORI_allB;
        margSF_allBB(iB_start:iB_end, :, :, :) = margSF_allB;
        % kernel fitting (ORI)
        margPred_ORI_allBB(iB_start:iB_end, :, :, :, :) = margPred_ORI_allB;        % 4D, ntypes x ncomb6 x nlines x nfilters_itp
        margParams_ORI_allBB(iB_start:iB_end, :, :, :, :) = margParams_ORI_allB; % 4D, ntypes x ncomb6 x nlines x nparams
        margR2_ORI_allBB(iB_start:iB_end, :, :, :) = margR2_ORI_allB;                 % 3D, ntypes x ncomb6 x nlines
        % kernel fitting (SF)
        margPred_SF_allBB(iB_start:iB_end, :, :, :, :) = margPred_SF_allB;
        margParams_SF_allBB(iB_start:iB_end, :, :, :, :) = margParams_SF_allB;
        margR2_SF_allBB(iB_start:iB_end, :, :, :) = margR2_SF_allB;
    case 6
        %% create containers for IDVD data, all batches [_allBB]
        % % linear regression slopes - behavioral
        pC_tgt_allBB = nan(nB, ntypes, nLoc8, nbins_e);
        ebin_tgt_allBB = pC_tgt_allBB;
        pC_tgt_allBB_norm = pC_tgt_allBB;
        ebin_tgt_allBB_norm = pC_tgt_allBB;
        
        % kernel
        kernels2D_allBB = nan(nB, ntypes, nLoc8, nfiltersOri, nfiltersSF); % ntypes, nLoc8, nfiltersOri, nfiltersSF
        intercept2D_allBB = kernels2D_allBB;
        margORI_allBB = nan(nB, ntypes, nfiltersOri, nLoc8);
        margSF_allBB = nan(nB, ntypes, nfiltersSF, nLoc8);
        R2_2D_allBB = nan(nB, ntypes, nLoc8, nfiltersOri, nfiltersSF);
        R2_Tjur_allBB = R2_2D_allBB;
        pValues_allBB = nan(nB, ntypes, nLoc8, nfiltersOri, nfiltersSF, 2); % ntypes, nLoc8, nfiltersOri, nfiltersSF, 2; % 2 = p of intercept+p of slope
        pCat_allBB = R2_2D_allBB;
        sep_allBB = nan(nB, ntypes, nLoc8); % ntypes, nLoc8
        
        % ORI-kernel fitting
        margPred_ORI_allBB = nan(nB, ntypes, ncomb6, nlines, nfiltersOri*nitp);
        margParams_ORI_allBB = nan(nB, ntypes, ncomb6, nlines, nparamsOri);
        margR2_ORI_allBB = nan(nB, ntypes, ncomb6, nlines);
        
        % SF - kernel fitting
        margPred_SF_allBB = nan(nB, ntypes, ncomb6, nlines, nfiltersSF*nitp);
        margParams_SF_allBB = nan(nB, ntypes, ncomb6, nlines, nparamsSF);
        margR2_SF_allBB = margR2_ORI_allBB;
    case 7
        %% save all subj
        % reorder the location (to select the dominant eye side for HM)
        if ~isnan(eyeDom)
            % slope
            pC_tgt_allBB = reorderLoc(pC_tgt_allBB, 1, eyeDom(isubj));
            ebin_tgt_allBB = reorderLoc(ebin_tgt_allBB, 1, eyeDom(isubj));
            pC_tgt_norm_allBB = reorderLoc(pC_tgt_norm_allBB, 1, eyeDom(isubj));
            ebin_tgt_norm_allBB = reorderLoc(ebin_tgt_norm_allBB, 1, eyeDom(isubj));
            % kernel
            kernels2D_allBB = reorderLoc(kernels2D_allBB, 2, eyeDom(isubj));
            intercept2D_allBB = reorderLoc(intercept2D_allBB, 2, eyeDom(isubj));
            R2_2D_allBB = reorderLoc(R2_2D_allBB, 2, eyeDom(isubj));
            R2_Tjur_allBB = reorderLoc(R2_Tjur_allBB, 2, eyeDom(isubj));
            pValues_allBB = reorderLoc(pValues_allBB, 3, eyeDom(isubj));
            pCat_allBB = reorderLoc(pCat_allBB, 2, eyeDom(isubj));
            sep_allBB = reorderLoc(sep_allBB, 4, eyeDom(isubj));
            % MARG
            margORI_allBB = reorderLoc(margORI_allBB, 5, eyeDom(isubj));
            margSF_allBB = reorderLoc(margSF_allBB, 5, eyeDom(isubj));
            % kernel fitting (ORI)
            margPred_ORI_allBB = reorderLoc(margPred_ORI_allBB, 6, eyeDom(isubj));
            margParams_ORI_allBB = reorderLoc(margParams_ORI_allBB, 6, eyeDom(isubj));
            margR2_ORI_allBB = reorderLoc(margR2_ORI_allBB, 7, eyeDom(isubj));
            % kernel fitting (SF)
            margPred_SF_allBB = reorderLoc(margPred_SF_allBB, 6, eyeDom(isubj));
            margParams_SF_allBB = reorderLoc(margParams_SF_allBB, 6, eyeDom(isubj));
            margR2_SF_allBB = reorderLoc(margR2_SF_allBB, 7, eyeDom(isubj));
        end
        
        % behav
        data_aveBlk_allSubj(isubj, :, :) = data_aveBlk_perComb;
        data_semBlk_allSubj(isubj, :, :) = data_semBlk_perComb;
        dprime_allT_allSubj(isubj, :) = dprime_allT_perComb;
        criterion_allT_allSubj(isubj, :) = criterion_allT_perComb;
        pA_allT_allSubj(isubj, :, :) = pA_allT_perComb;
        pC_allT_allSubj(isubj, :, :) = pC_allT_perComb;
        
        % slope
        pC_tgt_allSubj(isubj, :, :, :, :) = pC_tgt_allBB;
        ebin_tgt_allSubj(isubj, :, :, :, :) = ebin_tgt_allBB;
        pC_tgt_norm_allSubj(isubj, :, :, :, :) = pC_tgt_norm_allBB;
        ebin_tgt_norm_allSubj(isubj, :, :, :, :) = ebin_tgt_norm_allBB;
        % kernel
        kernels2D_allSubj(isubj, :, :, :, :, :) = kernels2D_allBB;
        intercept2D_allSubj(isubj, :, :, :, :, :) = intercept2D_allBB;
        R2_2D_allSubj(isubj, :, :, :, :, :) = R2_2D_allBB;
        R2_Tjur_allSubj(isubj, :, :, :, :, :) = R2_Tjur_allBB;
        pValues_allSubj(isubj, :, :, :, :, :, :) = pValues_allBB;
        pCat_allSubj(isubj, :, :, :, :, :) = pCat_allBB;
        sep_allSubj(isubj, :, :, :) = sep_allBB;
        % MARG
        margORI_allSubj(isubj, :, :, :, :) = margORI_allBB;
        margSF_allSubj(isubj, :, :, :, :) = margSF_allBB;
        % ORI - kernel fitting
        margPred_ORI_allSubj(isubj, :, :, :, :, :) = margPred_ORI_allBB;
        margParams_ORI_allSubj(isubj, :, :, :, :, :) = margParams_ORI_allBB;
        margR2_ORI_allSubj(isubj, :, :, :, :) = margR2_ORI_allBB;
        % SF - kernel fitting
        margPred_SF_allSubj(isubj, :, :, :, :, :) = margPred_SF_allBB;
        margParams_SF_allSubj(isubj, :, :, :, :, :) = margParams_SF_allBB;
        margR2_SF_allSubj(isubj, :, :, :, :) = margR2_SF_allBB;
end

%%
function data = reorderLoc(data, ID, eyeDom)

% eyeDom = 1: RIGHT Eye Dominant
% eyeDom = 0: LeftEye Dominant
% 2: Left [4,1] of combInd
% 4; Right [4,2]
% 6: HM ([2,1])

switch ID
    case 1 % pC/ebin_tgt(_norm)_perComb [3D] ntypes x nLoc8 x nbins
        if eyeDom, data(:, :, 6, :) = data(:, :, 4, :); else, data(:, :, 6, :) = data(:, :, 2, :); end
    case 2 % kernels2D_perComb, intercept2D_perComb, R2_2D_perComb, R2_Tjur_perComb, pCat_perComb
        % [4D] ntypes x nLoc8 x nfiltersOri x nfiltersSF
        if eyeDom, data(:, :, 6, :, :) = data(:, :, 4, :, :); else, data(:, :, 6, :, :) = data(:, :, 2, :, :); end
    case 3 % pValues_perComb [5D]
        if eyeDom, data(:, :, 6, :, :, :) = data(:, :, 4, :, :, :); else, data(:, :, 6, :, :, :) = data(:, :, 2, :, :, :); end
    case 4 % sep_perComb [2D] ntyps x nLoc8
        if eyeDom,
            data(:, :, 6) = data(:, :, 4);
        else, data(:, :, 6) = data(:, :, 2);
        end
    case 5 % margORI/SF_perComb [3D] ntypes x nfilters x nLoc8
        if eyeDom, data(:, :, :, 6) = data(:, :, :, 4); else, data(:, :, :, 6) = data(:, :, :, 2); end
    case 6 % margPred/Params_ORI/SF_perComb [4D] ntypes x ncomb6 x nlines x nfilters_itp/nparams
        if eyeDom, data(:, :, 2 ,1, :) = data(:, :, 4, 2, :); else, data(:, :, 2, 1, :) = data(:, :, 4, 1, :); end
    case 7 % margR2_ORI/SF_perComb, [3D] ntypes x ncomb6 x nlines
        if eyeDom, data(:, :, 2 ,1) = data(:, :, 4, 2); else, data(:, :, 2, 1) = data(:, :, 4, 1); end
end % switch ID

end


