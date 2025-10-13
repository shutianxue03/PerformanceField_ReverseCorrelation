
% generate 2D kernel matrix for each of the 4 categories: Hit, FA, Miss, CR
% then compute the Classification KERNELs

namesCategories = {'Hit', 'FA', 'Miss', 'CR'};
nc = length(namesCategories);
% itgtPrs, answer, correctness, RT, tgt phase

kernels2D_perLoc = nan(nLoc, nc, nfiltersOri, nfiltersSF);
R2_2D_kernel_perLoc = kernels2D_perLoc;
R2_Tjur_perLoc = R2_2D_kernel_perLoc ;
pValues_perLoc = nan(nLoc, nc, nfiltersOri, nfiltersSF, 2);
seprblity_kernel_perLoc = nan(nLoc, nc);
CK_perLoc = nan(nLoc, nfiltersOri, nfiltersSF); % classfication kernel

for iLoc = 1:nLoc
    for ic = 1:4
        switch ic
            case 1, indC = (data_both_B(:, iLoc, 1) == 1) & (data_both_B(:, iLoc, 2) == 1);
            case 2, indC = (data_both_B(:, iLoc, 1) == 0) & (data_both_B(:, iLoc, 2) == 1);
            case 3, indC = (data_both_B(:, iLoc, 1) == 1) & (data_both_B(:, iLoc, 2) == 0);
            case 4, indC = (data_both_B(:, iLoc, 1) == 0) & (data_both_B(:, iLoc, 2) == 0);
        end
        
        ntrials_c = sum(indC);
        e_iLoc = squeeze(energy2D_allT_perLoc(iLoc, indC, :, :));
        resp_iLoc = squeeze(data_both_B(indC, iLoc, 2));
        fprintf('r%s=%.2f\n', namesCategories{ic}, mean(resp_iLoc));
        
        [kernel2D_, ~, R2_2D_, R2_Tjur_, p_] = SX_sim07_RC(filtersSF_all, filtersOri_all, e_iLoc, resp_iLoc);
        kernels2D_perLoc(iLoc, ic, :, :) = kernel2D_;
        R2_2D_kernel_perLoc(iLoc, ic, :, :) = R2_2D_;
        R2_Tjur_perLoc(iLoc, ic, :, :) = R2_Tjur_;
        pValues_perLoc(iLoc, ic, :, :, :) = p_;
        seprblity_kernel_perLoc(iLoc, ic) = getSeparability(kernel2D_);
        
        fprintf('%s done\n', namesCategories{iLoc})
    end
    
    CK_perLoc(iLoc, :, :) = kernels2D_perLoc(iLoc, 1, :, :) + kernels2D_perLoc(iLoc, 2, :, :) - (kernels2D_perLoc(iLoc, 3, :, :)+kernels2D_perLoc(iLoc, 4, :, :));
    disp('%s done\')
end

%%
figure
for iLoc = 1:nLoc
    subplot(1,nLoc, iLoc)
    imagesc(squeeze(CK_perLoc(iLoc, :, :)))
end