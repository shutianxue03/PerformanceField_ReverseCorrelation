
xtext_p = [.58, 35]; ytext_p = [.015, .07]; % x and y of the text printing params

%%

kernels2D_allSubj_HM = kernels2D_allSubj_comp{6};
kernels2D_allSubj_VM = kernels2D_allSubj_comp{7};
kernels2D_allSubj_peri = kernels2D_allSubj_comp{8};

params_allLoc_PF3 = cell(1,2);


for ikernel = 2% 2 is ori, 1 is SF
    
    switch ikernel
        case 1
            xaxis = filtersSF_all_log_cut;
            nfilters = nfiltersSF_cut;
            idim = 3;
            modelInd = modelSF;
            nparams = nparamsSF;
            modelParamsNames = modelParamsNamesSF;
        case 2
            xaxis = filtersOri_all-90;
            nfilters = nfiltersOri;
            idim = 4;
            modelInd = modelOri;
            nparams = nparamsOri;
            modelParamsNames = modelParamsNamesOri;
    end
    
    
    % get mean and error
    if nsubj>1
        kernels_m_HM = squeeze(mean(kernels2D_allSubj_HM, idim));
        kernels_m_VM = squeeze(mean(kernels2D_allSubj_VM, idim));
        kernels_m_peri = squeeze(mean(kernels2D_allSubj_peri, idim));
        
        % mean
        kernels_ave_HM = squeeze(mean(kernels_m_HM, 1));
        kernels_ave_VM = squeeze(mean(kernels_m_VM, 1));
        kernels_ave_peri = squeeze(mean(kernels_m_peri, 1));
        
        % error
        kernels_err_HM = squeeze(std(kernels_m_HM, 1))/sqrt(nsubj);
        kernels_err_VM = squeeze(std(kernels_m_VM, 1))/sqrt(nsubj);
        kernels_err_peri = squeeze(std(kernels_m_peri, 1))/sqrt(nsubj);
    else
        kernels_m_HM = kernels2D_allSubj_HM;
        kernels_m_VM = kernels2D_allSubj_VM;
        kernels_m_peri = kernels2D_allSubj_peri;
        
        kernels_ave_HM = squeeze(mean(kernels2D_allSubj_HM, idim-1));
        kernels_ave_VM = squeeze(mean(kernels2D_allSubj_VM, idim-1));
        kernels_ave_peri = squeeze(mean(kernels2D_allSubj_peri, idim-1));
        
        kernels_err_HM = zeros(ntypes, nfilters);
        kernels_err_VM = kernels_err_HM;
        kernels_err_peri = kernels_err_HM;
    end
    
    params_allLoc_ = nan(length(PFnames), nparams);
    
     if publishFlag, figure('Position', [0 0 length(PFnames)*400 300]), end
     
    for ititle = 2% : length(PFnames)
        switch ititle
            case 1, kernels_ave = kernels_ave_HM; kernels_err = kernels_err_HM;
            case 2, kernels_ave = kernels_ave_VM; kernels_err = kernels_err_VM;
            case 3, kernels_ave = kernels_ave_peri; kernels_err = kernels_err_peri;
        end
        
        if publishFlag, subplot(1,length(PFnames),ititle), hold on
        else, figure('Position', [0 0 400 300]), hold on
        end
        
        for itype = 3%1:ntypes
            % data
            plot(xaxis, kernels_ave(itype, :), 'o', 'color', colorsType{itype})
            
            % error
            if nsubj>1
                ub = kernels_ave(itype, :) +  kernels_err(itype, :);
                lb = kernels_ave(itype, :) -  kernels_err(itype, :);
                x_flip = [xaxis, flip(xaxis)];
                y = [lb, flip(ub)];
                patch('Faces',1:2*nfilters,'Vertices', [x_flip', y'], 'FaceColor', colorsType{itype}, 'FaceAlpha', faceAlpha,'EdgeColor','none');
            end
            
            % fitting
            if ikernel==1, [kernels_pred, params, Rsquared, ~] = SX_sim08_fit(2.^xaxis, kernels_ave(itype, :), modelSF, fitMode);
            else, [kernels_pred, params, Rsquared, ~] = SX_sim08_fit(xaxis, kernels_ave(itype, :), modelOri, fitMode);
            end
            params_allLoc_(ititle, :) = params;
            % plot fitting
            plot(interp(xaxis, nx_itpl), interp(kernels_pred, nx_itpl), 'r-')
            % print text
            params_text = cell(1, nparams);
            for iparam = 1:nparams, params_text{iparam} = sprintf('%s = %.2f', modelParamsNames{iparam}, params(iparam)); end
%             if ikernel == 1, annotation('textbox', [0.3, 0.3, 0.1, 0.1], 'String', params_text,'FitBoxToText','on', 'fontsize', 12, 'BackgroundColor', 'w')
%             else, annotation('textbox', [0.65, 0.55, 0.1, 0.1], 'String', params_text,'FitBoxToText','on', 'fontsize', 12, 'BackgroundColor', 'w')
%             end
            
            % print R2
            annotation('textbox', [0.175, .8, 0.1, 0.1], 'String',  sprintf('R^2 = %.3f', Rsquared),'fontsize', sz_text, 'color', 'r', 'EdgeColor', 'w')
            
        end
        
        diffFlag = 1; RCplot_marg_format
         
        % save fig
        if ~publishFlag      
        if nsubj==1, saveas(gcf, sprintf('publishedPDFs/fig/%s_%s_marg_%s.jpg', subjName, xlabels{ikernel}(1:2), PFnames{ititle}))
        else, saveas(gcf, sprintf('publishedPDFs/fig/n%d_%s_marg_%s.jpg', nsubj, xlabels{ikernel}(1:2), PFnames{ititle}))
        end
        end
        
        params_allLoc_PF3{ikernel} = params_allLoc_;
    end
    
    % save fitting params
    if nsubj>1, save(sprintf('Data/params_n%d_PF%d',nsubj, length(PFnames)), 'params_allLoc_PF3')
    else, save(sprintf('Data/params_%s_PF%d', subjName,length(PFnames)), 'params_allLoc_PF3')
    end
    
    %% plot diff between marginalized kernels
    %     kernels_comp = cell(length(diffTypeNames), 1);
    %
    %     for iLoc = 1:nLoc, kernels_comp{iLoc} = squeeze(kernels_ave(:,iLoc, :)); end
    %     kernels_comp{6} = kernels_m_HM;
    %     kernels_comp{7} = kernels_m_VM;
    %     kernels_comp{8} = kernels_m_peri;
    %
    %     kernels_comp_allSubj = cell(length(diffTypeNames), 1);
    %     for iLoc = 1:nLoc, kernels_comp_allSubj{iLoc} = squeeze(kernels_allSubj(:,:, iLoc, :)); end
    %     kernels_comp_allSubj{6} = squeeze(mean(kernels_allSubj(:,:,[2,4],:),3));
    %     kernels_comp_allSubj{7} = squeeze(mean(kernels_allSubj(:,:,[3,5],:),3));
    %     kernels_comp_allSubj{8} = squeeze(mean(kernels_allSubj(:,:,[2:5],:),3));
    %
    %     %%
    %     figure('Position', [0 200 900 300])
    %     for icomp = 1:nComp
    %         diffType1 = diffType1_all(icomp);
    %         diffType2 = diffType2_all(icomp);
    %         subplot(1,3, icomp), hold on
    %
    %         for itype = 3%:ntypes
    %             kernels_m1 = kernels_comp{diffType1}(itype, :);
    %             kernels_m2 = kernels_comp{diffType2}(itype, :);
    %             kernels_m1_allSubj = squeeze(kernels_comp_allSubj{diffType1}(:, itype, :));
    %             kernels_m2_allSubj = squeeze(kernels_comp_allSubj{diffType2}(:, itype, :));
    %
    %             if nsubj>1, kernels_diff_err = std(kernels_m1_allSubj - kernels_m2_allSubj)/sqrt(nsubj);
    %             else, kernels_diff_err = zeros(1, nfilters);
    %             end
    %
    %             errorbar(xaxis, kernels_m1-kernels_m2, kernels_diff_err, 'o-', 'color', colorsType{itype})
    %             axis square
    %             if itype == 1, title(sprintf('%s vs. %s', diffTypeNames{diffType1}, diffTypeNames{diffType2})), end
    %         end
    %         yline(0);
    %         xline(x_xlline(ikernel)); % draw a vertical line
    %         xlabel(xlabels{ikernel})
    %         xticks(xticks_{ikernel})
    %         xticklabels(xticklabels_{ikernel})
    %         ylabel('kernel (a.u.)')
    %         ylim([yrange_diff{ikernel}(1), yrange_diff{ikernel}(end)])
    %         yticks(yrange_diff{ikernel})
    %         title(sprintf('%s vs. %s', diffTypeNames{diffType1}, diffTypeNames{diffType2}))
    %     end
    %
    %     set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)
    
end

