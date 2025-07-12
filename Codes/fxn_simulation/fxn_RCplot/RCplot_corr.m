
% correlational analysis
% assess the correlation between the contrast thresh vs.. each parameter
icomb5 = [1,4,3,4,3];
iline5 = [1,1,2,2,1];
if nsubj>1, zeroMean_all = [0,1];
else, zeroMean_all = 0;
end

% get data index to conduct partial corr (controlling loc to reveal observer effect)
pcorr_indLoc = repmat(1:nLoc5, nsubj,1);
pcorr_indLoc = pcorr_indLoc(:);


%%
for zeroMean = zeroMean_all
    for ix = 1:2 % ix = 1: CST in on xaxis; ix=2: pA is on xaxis
        % extract x (CST or pA)
        switch ix
            case 1
                xlabel_ = 'CST thresh';
            case 2
                xlabel_ = 'Response consistency (pA)';
        end
        % get behav data 
        bm_aveBlk_perLoc = squeeze(data_aveBlk_allSubj(:, ix+3, 1:nLoc5));
        bm_semBlk_perLoc = squeeze(data_semBlk_allSubj(:, ix+3, 1:nLoc5));
        
        for ifeature = 1:2
            ifamily = ifamily_perF(ifeature);
            nameParams = namesParams_all{ifamily};
            nparams = length(nameParams);
            
            figure('Position', [zeroMean*(nparams*300) 0 nparams*300 1000])
            
            for itype = 1:ntypes
                
                for iparam = 1:nparams
                    
                    % reorganize params (4x2 to 5x1)
                    params_med_perLoc = nan(nsubj, nLoc5);
                    params_CI_neg_perLoc = params_med_perLoc;
                    params_CI_pos_perLoc = params_med_perLoc;
                    for iLoc = 1:nLoc5
                        icomb = combInd_inv(iLoc, 1);
                        iline = combInd_inv(iLoc, 2);
                        if nsubj==1
                            params_med_perLoc(:, iLoc) = params_medB_allSubj_all{ifeature, itype, iparam, icomb}(iline);
                            params_CI_neg_perLoc(:, iLoc) = params_CI_B_neg_allSubj_all{ifeature, itype, iparam, icomb}(iline);
                            params_CI_pos_perLoc(:, iLoc) = params_CI_B_pos_allSubj_all{ifeature, itype, iparam, icomb}(iline);
                        else
                            params_med_perLoc(:, iLoc) = params_medB_allSubj_all{ifeature, itype, iparam, icomb}(:, iline);
                            params_CI_neg_perLoc(:, iLoc) = params_CI_B_neg_allSubj_all{ifeature, itype, iparam, icomb}(:, iline);
                            params_CI_pos_perLoc(:, iLoc) = params_CI_B_pos_allSubj_all{ifeature, itype, iparam, icomb}(:, iline);
                        end
                    end
                    
                    % zero-mean params
                    if zeroMean
                        params_med_perLoc = params_med_perLoc - mean(params_med_perLoc, 1);
                        bm_aveBlk_perLoc = bm_aveBlk_perLoc - mean(bm_aveBlk_perLoc, 1);
                    end
                    
                    subplot(ntypes, nparams, (itype-1)*nparams+iparam)
                    hold on, box on, grid on
                    text_corr = cell(1,nLoc5);
                    
                    for iLoc = 1:nLoc5
                        % plot IDVD data
                        if nsubj>1
                        errorbar(bm_aveBlk_perLoc(:, iLoc), params_med_perLoc(:,iLoc), ...
                            params_CI_neg_perLoc(:, iLoc), params_CI_pos_perLoc(:, iLoc), ...
                            bm_semBlk_perLoc(:, iLoc), bm_semBlk_perLoc(:, iLoc), ...
                            'o', 'color', colors5Loc(iLoc, :), 'CapSize', 0)
                        else
                            errorbar(bm_aveBlk_perLoc(iLoc), params_med_perLoc(iLoc), ...
                            params_CI_neg_perLoc(iLoc), params_CI_pos_perLoc(iLoc), ...
                            bm_semBlk_perLoc(iLoc), bm_semBlk_perLoc(iLoc), ...
                            'o', 'color', colors5Loc(iLoc, :), 'CapSize', 0)
                        end
                    end % end of iLoc
                    
                    % linear regression on data of ALL loc
                    title_corr = []; % do not delete
                    if nsubj > 1
                        % Pearson's corr
                        x = bm_aveBlk_perLoc(:);
                        y = params_med_perLoc(:);
                        
                        [r, p] = corr(x, y);
                        lm = polyfit(x, y, 1);
                        yfit = polyval(lm, x);
                        plot(x, yfit, 'k', 'handlevisibility', 'off', 'linewidth', 2)
                        R2 = getR2(y, yfit);
                        ss = getString_starts(p);
                        
                        % partial corr (controlling for location)
                        if zeroMean
                            r_loc = nan; ss_loc = nan;
                        else
                            [r_loc, p_loc] = partialcorr([x, y], pcorr_indLoc);
                            r_loc = r_loc(2,1);
                            p_loc = p_loc(2,1);
                            ss_loc = getString_starts(p_loc);
                        end
                        
                        % make the figure title (corr coefficient and p values)
                        title_corr = sprintf('r=%.2f%s (r=%.2f%s)', r, ss, r_loc, ss_loc);
                        
                    end
                    
                    % extra lines
                    if ifeature == 1 && iparam == 3, yline(0, 'color', ones(1,3)*.5); end
                    
                    % labels and title
                    %                 if (iparam == nparams) && (itype == 1), xlabel(xlabel_), end
                    %                 if nsubj>1, if (iparam == 1) && (itype == 1), legend(subjList, 'NumColumns', 2, 'Location', 'southwest'); end, end
                    title(sprintf('%s\n%s', namesType{itype}, title_corr))
                    ylabel(sprintf('%s %s', namesFeature{ifeature}, nameParams{iparam}))
                    
                    axis square
                end % end of iparam
            end % end of itype
            
            set(findall(gcf, '-property', 'FontSize'), 'FontSize',18)
            make_sgtitle(sprintf('[%s] params vs. %s', namesFeature{ifeature}, xlabel_), subjName, nB, nsubj, nAllTrials)
            
        end % end of ifeature
    end % end of ix
end % end of zeroMean = [0,1]
