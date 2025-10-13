
%% Bootstrap the fitted line to the tuning

if nsubj==1, RCbootName = sprintf('Data/%s/%s%d_RCboot%s.mat', subjName, subjName, nblocks, title_);
else, RCbootName = sprintf('Data/RCboot_n%d_%s.mat', nsubj, title_);
end
RCbootDir = dir(RCbootName);

%%
if ~isempty(RCbootDir), load(RCbootName), fprintf('boot file exists.\n')
    
else
    fprintf('boot on model fitting starts (nB=%d).\n', nB_k)
    
    % empty containers
    kpred_allB_allLines = cell(ntypes_marg, 2, ncomp); % PRS/ABS/BOTH, ORI/SF, ncomp
    kparams_allB_allLines = kpred_allB_allLines;
    kR2_allB_allLines = kpred_allB_allLines;
    marg_ave_allLines = kpred_allB_allLines;
    marg_SEM_allLines = kpred_allB_allLines;
    
    for itype_used = 1:ntypes_marg
        
        for ifeature = 1:2
            xaxis = axis_all{ifeature};
            nfilters = length(xaxis);
            nameModelParams = namesModelParams{ifeature};
            nparams = length(nameModelParams);
            
            switch ifeature
                case 1, marg_comp_allSubj = dataOri_allSubj_comp; modelKernel = modelOri;
                case 2, marg_comp_allSubj = dataSF_allSubj_comp; modelKernel = modelSF;
            end
            
            for icomp = 1:ncomp
                kpred_allB = nan(2, nB_k, nfilters); % 1=line #1, 2 = line#2 (two lines to be compared)
                kparams_allB = nan(2, nB_k, nparams);
                kR2_allB = nan(2, nB_k);
                
                for iline = 1:2
                    % group average
                    if nsubj > 1 % weighted average
                        marg_comp = marg_comp_allSubj{compInd(icomp, iline)}(:, itype_used, :);
                        marg_ave = 0;
                        for isubj = 1:nsubj, marg_ave = marg_ave + squeeze(marg_comp(isubj, :, :)) * ntrialsProp(isubj); end
                        marg_SEM = squeeze(std(marg_comp, [], 1)/sqrt(nsubj));
                    else
                        marg_comp = marg_comp_allSubj{compInd(icomp, iline)}(itype_used, :);
                        marg_ave = squeeze(marg_comp);
                        marg_SEM = zeros(size(marg_ave));
                    end
                    
                    % to ensure the correct shape
                    marg_ave = reshape(marg_ave, 1, length(marg_ave));
                    marg_SEM = reshape(marg_SEM, 1, length(marg_SEM));
                    
                    % save 
                    marg_ave_allLines{itype_used, ifeature, icomp}(iline, :) = marg_ave;
                    marg_SEM_allLines{itype_used, ifeature, icomp}(iline, :) = marg_SEM;
                    
                    % boostrapping
                    parfor iB = 1:nB_k
                        dataInd = randi(nfilters, 1, nfilters); % resample the marginalized kernel
                        switch ifeature
                            case 1
                                [~, params_est, Rsquared, ~] = SX_sim08_fit(xaxis(dataInd), marg_ave(dataInd), modelKernel, fitMode);
                                kpred = predSFkernel(xaxis, modelKernel, params_est, 0);
                            case 2
                                [~, params_est, Rsquared, ~] = SX_sim08_fit(2.^xaxis(dataInd), marg_ave(dataInd), modelKernel, fitMode);
                                kpred = predSFkernel(2.^xaxis, modelKernel, params_est, 0);
                        end
                        kpred_allB(iline, iB, :) = kpred;
                        kparams_allB(iline, iB, :) = params_est;
                        kR2_allB(iline, iB) = Rsquared;
                    end
                    
                    % save to the bigger structure
                    kpred_allB_allLines{itype_used, ifeature, icomp} = kpred_allB;
                    kparams_allB_allLines{itype_used, ifeature, icomp} = kparams_allB;
                    kR2_allB_allLines{itype_used, ifeature, icomp} = kR2_allB;
                    
                end
            end
        end
    end
    
    fprintf('bootstrapping on fitting done\n')
    % save the bootstrapped data
    save(RCbootName, 'kpred_allB_allLines', 'kparams_allB_allLines', 'kR2_allB_allLines', 'marg_ave_allLines', 'marg_SEM_allLines')
end

