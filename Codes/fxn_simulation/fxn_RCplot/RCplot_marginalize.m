
clear patch

xtext_p = [.58, 35]; 
ytext_p = [0, .07]; % x and y of the text printing params for SF and ori

% marginalize for all subjects
kernelsSF_allSubj = squeeze(mean(kernels2D_allSubj_mirror, 4));
kernelsOri_allSubj = squeeze(mean(kernels2D_allSubj_mirror, 5));

%% get mean and error
RC_getKernelMeanError

params_allLoc_loc5 = cell(1,2);

%%
for ikernel = 2 % SF and ori
    switch ikernel
        case 1
            kernels_ave_allLoc = kernelsSF_ave;
            kernels_err_allLoc = kernelsSF_err;
            xaxis = filtersSF_all_log_cut;
            modelInd = modelSF;
            nparams = nparamsSF;
            modelParamsNames = modelParamsNamesSF;
        case 2
            kernels_ave_allLoc = kernelsOri_ave;
            kernels_err_allLoc = kernelsOri_err;
            xaxis = filtersOri_all-90;
            modelInd = modelOri;
            nparams = nparamsOri;
            modelParamsNames = modelParamsNamesOri;
    end
    nfilters = length(xaxis);
    params_allLoc_ = nan(nLoc, nparams);
    
    if publishFlag, figure('Position', [0 0 1000 800]), end
    
    for iLoc = 3%:nLoc
        % plot
        if publishFlag, subplot(3,3, subplot_locs(iLoc)), hold on
        else, figure('Position', [0 0 400 300]), hold on
        end
        
        for itype = 3%1:ntypes
            kernels_ave =  squeeze(kernels_ave_allLoc(itype, iLoc, :));
            kernels_err = squeeze(kernels_err_allLoc(itype, iLoc, :));
            
            % data
            plot(xaxis, kernels_ave.', 'o', 'color', colorsType{itype})
            
            % error
            if nsubj>1
                ub = kernels_ave + kernels_err;
                lb = kernels_ave - kernels_err;
                x = [xaxis, flip(xaxis)];
                y = [lb; flip(ub)];
                patch('Faces',1:2*nfilters,'Vertices', [x', y], 'FaceColor', colorsType{itype}, 'FaceAlpha', faceAlpha,'EdgeColor','none');
            end
            
            % fitting
            if ikernel==1, [kernels_pred, params, Rsquared, ~] = SX_sim08_fit(2.^xaxis, reshape(kernels_ave, 1, length(xaxis)), modelSF, fitMode);
            else, [kernels_pred, params, Rsquared, ~] = SX_sim08_fit(xaxis, reshape(kernels_ave, 1, length(xaxis)), modelOri, fitMode);
            end
            params_allLoc_(iLoc, :) = params;
            
            % plot fitting
            plot(interp(xaxis, nx_itpl), interp(kernels_pred, nx_itpl), '-r')
            
            % print text
            params_text = cell(1, nparams);
            for iparam = 1:nparams, params_text{iparam} = sprintf('%s = %.2f', modelParamsNames{iparam}, params(iparam));end
            % text of params
%             if ikernel == 1, annotation('textbox', [0.3, 0.3, 0.1, 0.1], 'String', params_text,'FitBoxToText','on', 'fontsize', 12, 'BackgroundColor', 'w')
%             else, annotation('textbox', [0.65, 0.55, 0.1, 0.1], 'String', params_text,'FitBoxToText','on', 'fontsize', 12, 'BackgroundColor', 'w')
%             end
            % print R2
            text(xaxis(1), yrange{ikernel}(end)-.01, sprintf('R^2 = %.3f', Rsquared), 'FontSize', sz_text, 'color', 'r') % text of R2
            
            diffFlag = 0; RCplot_marg_format
        end
        
        % save fig
        if ~publishFlag
            if nsubj==1, saveas(gcf, sprintf('publishedPDFs/fig/%s_%s_marg_%s.jpg', subjName, xlabels{ikernel}(1:2), locNames{iLoc}))
            else, saveas(gcf, sprintf('publishedPDFs/fig/n%d_%s_marg_%s.jpg', nsubj, xlabels{ikernel}(1:2), locNames{iLoc}))
            end
        end
    end
    params_allLoc_loc5{ikernel} = params_allLoc_;
end

% save fitting params
if nsubj>1, save(sprintf('Data/params_n%d_loc%d',nsubj, nLoc), 'params_allLoc_loc5')
else, save(sprintf('Data/params_%s_loc%d', subjName,nLoc), 'params_allLoc_loc5')
end

