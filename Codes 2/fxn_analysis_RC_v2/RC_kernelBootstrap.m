
% bootstrap kernels generated from PRS, ABS and BOTH 

fprintf('boostrapping kernels ...\n')
nB = 5; % number of bootstrap rounds
prop = 0.8; % prop of trials 
CI_level = .68; % confidence interval (plotted as errorbar)  (Fernandez, 2021)

ntrialsB = ntrials * prop;

kernelOri_boot = nan(ntypes, nLoc, nB, nfiltersOri);
kernelSF_boot = nan(ntypes, nLoc, nB, nfiltersSF);
kernelOri_CI68 = nan(2, ntypes, nLoc, nfiltersOri); % 2: upper and lower bound; 3: prs,abs and both
kernelOri_CI68_errbar = kernelOri_CI68; % 2: upper and lower bound; 3: prs,abs and both
kernelSF_CI68 = nan(2, ntypes, nLoc, nfiltersSF); % 2: upper and lower bound; 3: prs,abs and both
kernelSF_CI68_errbar = kernelSF_CI68; % 2: upper and lower bound; 3: prs,abs and both

for itype = 3 % prs, abs, both
    fprintf('\n%s trials ...', typeNames{itype})
    for iLoc = 1:nLoc
        for iB = 1:nB
            % get index
            switch itype
                case 1, ind = randi(ntrials, 1, ntrialsB);
                case 2, ind = randi(ntrials, 1, ntrialsB) + ntrials;
                case 3, ind = [randi(ntrials, 1, ntrialsB), randi(ntrials, 1, ntrialsB)+ntrials];
            end
            
            % compute kernel
            [kernel2D_boot, var_exp] = SX_sim07_RC(filtersSF_all, filtersOri_all, squeeze(energy_norm_allLoc{iLoc}(ind, :, :)), data_both(ind, iLoc, 2));
            
            % mirror
            kernel2D_mirror_boot = getMirror(kernel2D_boot, 0, 1);
            
            % marginalize
            kernelSF_boot_ = mean(kernel2D_mirror_boot, 1);
            kernelOri_boot_ = mean(kernel2D_mirror_boot, 2).';
            kernelSF_boot(itype, iLoc, iB, :) = kernelSF_boot_;
            kernelOri_boot(itype, iLoc, iB, :) = kernelOri_boot_;
        end
        
        % get 68% CI
        kk = squeeze(kernelSF_boot(itype,iLoc,:,:));
        CI = quantile(kk, [.5-CI_level/2, .5+CI_level/2]);
        kernelSF_CI68(:, itype, iLoc, :) = CI;
        kernelSF_CI68_errbar(:, itype, iLoc, :) = [mean(kk) - CI(1,:); CI(2,:) - mean(kk)];
        
        kk = squeeze(kernelOri_boot(itype,iLoc,:,:));
        CI = quantile(kk, [.5-CI_level/2, .5+CI_level/2]);
        kernelOri_CI68(:, itype, iLoc, :) = CI;
        kernelOri_CI68_errbar(:, itype, iLoc, :) = [mean(kk) - CI(1,:); CI(2,:) - mean(kk)];
        
        fprintf('L%d done. ', iLoc)
    end
end



