
% cut ORI and SF because kernel is messy
itype = 2;
% cut
kernels2D_allB = kernels2D_allB(:, :, :, cut_ORI, cut_SF);

%% mirror
if flag_mirrorMapping
    indMir = (nORI-1)/2;
    kk_left = kernels2D_allB(:, :, :, 1:indMir, :);
    kk_right = kernels2D_allB(:, :, :, indMir+2:end, :);
    kk_mid = kernels2D_allB(:, :, :, indMir+1,:);
    kk_ave = (kk_left + flip(kk_right, 4))/2;
    kernels2D_allB_mir = cat(4, kk_ave, kk_mid, flip(kk_ave, 4));
    
    % get similarity after mirroring
    for iB = 1:nB
        for iiLoc = 1:2
            similarityAfterMirroring_allB(iB, iiLoc) = corr2(squeeze(kernels2D_allB(iB, iiLoc, itype, :, :)), squeeze(kernels2D_allB_mir(iB, iiLoc, itype, :, :)));
            %             figure
            %             subplot(1,2,1), imagesc(squeeze(kernels2D_allB(iB, iiLoc, itype, :, :))), axis square, colorbar
            %             subplot(1,2,2), imagesc(squeeze(kernels2D_allB_mir(iB, iiLoc, itype, :, :))), axis square, colorbar
        end
    end
    
    for iiLoc = 1:2
        kk = squeeze(median(kernels2D_allB(:, iiLoc, itype, :, :)));
        kk_m = squeeze(median(kernels2D_allB_mir(:, iiLoc, itype, :, :)));
        if flag_plot
        figure
        subplot(2,2,1), imagesc(kk), axis square, colorbar, yline((nORI+1)/2, 'r-'); xline((nSF+1)/2, 'r-'); title('Raw')
        subplot(2,2,2), imagesc(kk_m), axis square, colorbar, yline((nORI+1)/2, 'r-'); xline((nSF+1)/2, 'r-'); title('Mirrored')
        subplot(2,2,3), hold on
        plot(axis_tuning{1}, mean(kk, 2), '-')
        plot(axis_tuning{1}, mean(kk_m, 2), '--')
        xlabel('ORI')
        subplot(2,2,4), hold on
        plot(axis_tuning{2}, mean(kk, 1), '-')
        plot(axis_tuning{2}, mean(kk_m, 1), '--')
        xlabel('SF')
        sgtitle(namesLocComb{iLocComb_all(iiLoc)})
        end
    end
    kernels2D_allB = kernels2D_allB_mir;
end
median(similarityAfterMirroring_allB, 1)

%%
for iB =1:nB
    margORI_perComb = nan(2, ntypes, nORI);
    margSF_perComb = nan(2, ntypes, nSF);
    for iiLoc = 1:2
        e2D_cut_ = squeeze(kernels2D_allB(iB, iiLoc, itype, :, : ));
        e_min = min(e2D_cut_(:));
        e2D_cut = e2D_cut_ - e_min + eps;
        margORI = mean(e2D_cut, 2);
        margSF = mean(e2D_cut, 1);
        e2D_recon = mtimes(margORI, margSF);
        %         figure
        %         subplot(2,3,1), imagesc(e2D_cut_), colorbar
        %         subplot(2,3,2), plot(mean(e2D_cut_, 2)), [r,p] = corr(mean(e2D_cut_, 2), mean(e2D_cut, 2)); title(sprintf('corr=%.2f, p=%.3f', r,p))
        %         subplot(2,3,3), plot(mean(e2D_cut_, 1)), [r,p] = corr(mean(e2D_cut_, 1)', mean(e2D_cut, 1)'); title(sprintf('corr=%.2f, p=%.3f', r,p))
        %         subplot(2,3,4), imagesc(e2D_cut), colorbar
        %         subplot(2,3,5), plot(mean(e2D_cut, 2))
        %         subplot(2,3,6), plot(mean(e2D_cut, 1))
        
        sep_allB(iB, iiLoc, itype) = corr2(e2D_cut, e2D_recon);
        
        % marg
        margORI_allB(iB, iiLoc, itype, :) = margORI + e_min;
        margSF_allB(iB, iiLoc, itype, :) = margSF + e_min;
        margORI_perComb(iiLoc, :, :) = repmat(margORI + e_min, [1, ntypes]).'; % stupid code, need to feed three types into SX_RC7_fitting
        margSF_perComb(iiLoc, :, :) = repmat( margSF + e_min, [ntypes, 1]);
    end % iiLoc
    
    %%     fit marg and R2
    iF = 1; [margORI, margPred_ORI, margParams_ORI, margR2_ORI] = SX_RC7_fitting(iF, ...
        margORI_perComb, ub_full_all{ifamily_perF(iF)}, lb_full_all{ifamily_perF(iF)}, ...
        ifamily_perF, paramInd_perF, flag_standEnergy);
    %     %%
    %     x_intp = linspace(axis_tuning{1}(1), axis_tuning{1}(end), flag_interpolate);
    %     figure
    %     iiLoc = 1; subplot(1,2,1), hold on, plot(axis_tuning{1}, squeeze(margORI_perComb(iiLoc, 2, :)), 'o'), plot(x_intp, squeeze(margPred_ORI(iiLoc, 2, :)), '-'), xline(0, 'r-');
    %     iiLoc = 2; subplot(1,2,2), hold on, plot(axis_tuning{1}, squeeze(margORI_perComb(iiLoc, 2, :)), 'o'), plot(x_intp, squeeze(margPred_ORI(iiLoc, 2, :)), '-'), xline(0, 'r-');
    %         iiLoc = 1; subplot(1,2,1), hold on, plot(x_intp, squeeze(margORI(iiLoc, 2, :)), '.'), plot(x_intp, squeeze(margPred_ORI(iiLoc, 2, :)), '-'), xline(0, 'r-');
    %         iiLoc = 2; subplot(1,2,2), hold on, plot(x_intp, squeeze(margORI(iiLoc, 2, :)), '.'), plot(x_intp, squeeze(margPred_ORI(iiLoc, 2, :)), '-'), xline(0, 'r-');
    %%
    iF = 2; [margSF, margPred_SF, margParams_SF, margR2_SF] = SX_RC7_fitting(iF, ...
        margSF_perComb, ub_full_all{ifamily_perF(iF)}, lb_full_all{ifamily_perF(iF)}, ...
        ifamily_perF, paramInd_perF, flag_standEnergy);
    
    %%
    %     x_intp = linspace(axis_tuning{2}(1), axis_tuning{2}(end), flag_interpolate);
    %     figure
    %     iiLoc = 1; subplot(1,2,1), hold on, plot(axis_tuning{2}, squeeze(margSF_perComb(iiLoc, 2, :)), 'o'), plot(x_intp, squeeze(margPred_SF(iiLoc, 2, :)), '-'), xline(0, 'r-');
    %     xline(log2(margParams_SF(iiLoc, 2, 1)));
    %     xline(log2(margParams_SF(iiLoc, 2, 5)));
    %     iiLoc = 2; subplot(1,2,2), hold on, plot(axis_tuning{2}, squeeze(margSF_perComb(iiLoc, 2, :)), 'o'), plot(x_intp, squeeze(margPred_SF(iiLoc, 2, :)), '-'), xline(0, 'r-');
    %         xline(log2(margParams_SF(iiLoc, 2, 1)));
    %     xline(log2(margParams_SF(iiLoc, 2, 5)));
    %%
    margPredORI_allB(iB, :, :, :) = margPred_ORI;
    margPredSF_allB(iB, :, :, :) = margPred_SF;
    
    margR2ORI_allB(iB, :, :) = margR2_ORI;
    margR2SF_allB(iB, :, :) = margR2_SF;
    % get tuning
    for iiLoc = 1:2 % 2 locations
        for itype = 2%1:ntypes
            ifeature= 1;
            margTuningC_ORI = fxn_getTuningC(axis_tuning{ifeature}, ifeature, ifamily_perF(ifeature), squeeze(margPred_ORI(iiLoc, itype,:)), squeeze(margParams_ORI(iiLoc, itype,:)));
            ifeature= 2;
            margTuningC_SF = fxn_getTuningC(axis_tuning{ifeature}, ifeature, ifamily_perF(ifeature), squeeze(margPred_SF(iiLoc, itype,:)), squeeze(margParams_SF(iiLoc, itype,:)));
            % compile
            margTuningC_ORI_allB(iB, iiLoc, itype, :) = margTuningC_ORI;
            margTuningC_SF_allB(iB, iiLoc, itype, :) = margTuningC_SF;
        end % itype
    end % ii
    fprintf('.')
    
    
end % iBs

