
% close all
clear ydiffErr

%% settings
wd_border = 1;
sz_marker = 100;

% Bar-ORI
yticks_bar{1} = [0.02:.05:.22; 4:4:20; -.04:.015:.02];
ylim_bar{1} = [0, .24; 3.2, 20.8; -.046, .026];
axisTicks_scatter{1} = [.02:.1:.22; 8:7:22; -.04:.04:.04];
axisLim_scatter{1} = [0, .24; 6.6, 23.4; -.048, .048];
% Bar-SF
yticks_bar{2} = [1:.25:2; .01:.02:.09; 0:.1:.4; -.001:.0005:.001];
ylim_bar{2} = [.95, 2.05; 0 .1; 0, .4; -.001,.001];
axisTicks_scatter{2} = [1.5:.3:2.1; .01:.05:.11; 0:.2:.4;-.001:.001:.001];
axisLim_scatter{2} = [1.44, 2.16; 0, .12; 0, .4; -.001,.001];

%%
for ifeature = 1:2
    xaxis = axis_tuning{ifeature};
    nfilters = length(xaxis);
    xaxis_itp = axis_itp_tuning{ifeature};
    ifamily = ifamily_perF(ifeature);
    nameParams = namesParams_all{ifamily};
    nparams_full = length(nameParams);
    namesParams_unit = namesParams_unit_perF{ifeature};
    
    for itype = 3
        for iparam = 1:nparams_full
            
            figure('Position', [0 200 ncomb6*300 500]);
            
            for icomb = 1:ncomb6
                
                if paramInd_perF{ifeature}(iparam) == 1, nlines = 2; else, nlines = 1; end
                
                params_medB_allSubj = nan(nsubj, nlines);
                params_CI_B_neg_allSubj = params_medB_allSubj;
                params_CI_B_pos_allSubj = params_medB_allSubj;
                params_aveSubj  = nan(1, nlines);
                params_semSubj = params_aveSubj;
                
                for iline = 1:nlines
                    if (icomb == 2) && (iline == 1) && (eyeDom(1) > 0)
                        % right eye dominant, 6 (peri) is changed to 4 (RHM)
                        if eyeDom(isubj), icomb_real = 4; iline_real = 2;
                        else, icomb_real = 4; iline_real = 1;
                        end
                    else, icomb_real = icomb; iline_real = iline;
                    end
                    
                    if nsubj>1
                        params_medB_allSubj(:, iline) = params_medB_allSubj_all{ifeature, itype, iparam, icomb_real}(:, iline_real);
                        params_CI_B_neg_allSubj(:, iline) = params_CI_B_neg_allSubj_all{ifeature, itype, iparam, icomb_real}(:, iline_real);
                        params_CI_B_pos_allSubj(:, iline) = params_CI_B_pos_allSubj_all{ifeature, itype, iparam, icomb_real}(:, iline_real);
                    else
                        params_medB_allSubj(:, iline) = params_medB_allSubj_all{ifeature, itype, iparam, icomb_real}(iline_real, :);
                        params_CI_B_neg_allSubj(:, iline) = params_CI_B_neg_allSubj_all{ifeature, itype, iparam, icomb_real}(iline_real, :);
                        params_CI_B_pos_allSubj(:, iline) = params_CI_B_pos_allSubj_all{ifeature, itype, iparam, icomb_real}(iline_real, :);
                    end
                    params_aveSubj(iline) = params_aveSubj_all{ifeature, itype, iparam, icomb_real}(iline_real);
                    params_semSubj(iline) = params_semSubj_all{ifeature, itype, iparam, icomb_real}(iline_real);
                    
                    %%%%%%%%
                    %       bar       %
                    %%%%%%%%
                    subplot(2, ncomb6, icomb), hold on, box on
                    % group average
                    color = colors_comb(combInd(icomb_real, iline_real), :);
                    bar(iline, params_aveSubj(iline), 'FaceColor', 'w', 'EdgeColor', color, 'BarWidth', .5, 'linewidth', wd_border)
                    if nsubj == 1, errorbar(iline, params_medB_allSubj(iline), params_CI_B_neg_allSubj(iline), params_CI_B_pos_allSubj(iline), 'color', color, 'CapSize', 0);
                    else,errorbar(iline, params_aveSubj(iline), params_semSubj(iline), 'color', color, 'CapSize', 0)
                    end
                    
                    xticks(1:nlines)
                    xlim([.3, nlines+.7])
                    xticklabels({namesLocComb{combInd(icomb_real, 1:nlines)}}) % do NOT delete the extra curly brackets!!
                    
%                     yticks(yticks_bar{ifeature}(iparam, :))
%                     yticklabels(yticks_bar{ifeature}(iparam, :))
%                     ylim(ylim_bar{ifeature}(iparam, :))
                    
                end % end of iline
                
                %%%%%%%%
                %    scatter   %
                %%%%%%%%
                if nlines==2
                subplot(2, ncomb6, icomb+ncomb6), hold on, box on
                
                % diagonal line
                plot(axisLim_scatter{ifeature}, axisLim_scatter{ifeature}, 'k-', 'linewidth', wd_border)
                
                % IDVD data
%                 if nsubj>1
                    errorbar(params_medB_allSubj(:,1), params_medB_allSubj(:,2), params_CI_B_neg_allSubj(:,1), params_CI_B_pos_allSubj(:,1), params_CI_B_neg_allSubj(:,2), params_CI_B_pos_allSubj(:, 2), '.k', 'CapSize', 0, 'linewidth', wd_border/2)
                    scatter(params_medB_allSubj(:,1), params_medB_allSubj(:,2), sz_marker, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'linewidth', 2)
%                 else
%                     errorbar(params_medB_allSubj, params_medB_allSubj(:,2), params_CI_B_neg_allSubj(:,1), params_CI_B_pos_allSubj(:,1), params_CI_B_neg_allSubj(:,2), params_CI_B_pos_allSubj(:, 2), '.k', 'CapSize', 0, 'linewidth', wd_border/2)
%                     scatter(params_medB_allSubj, params_medB_allSubj(:,2), sz_marker, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'linewidth', 2)
%                 end
                % group ave (purple cross)
                if nsubj>1
                    errorbar(params_aveSubj(1), params_aveSubj(2), params_semSubj(1), params_semSubj(1), params_semSubj(2), params_semSubj(2), 'color', [1,0,1], 'linewidth', 2, 'CapSize', 0)
                end
                %
%                 xticks(axisTicks_scatter{ifeature}(iparam, :))
%                 yticks(axisTicks_scatter{ifeature}(iparam, :))
%                 
%                 xlim(axisLim_scatter{ifeature}(iparam, :))
%                 ylim(axisLim_scatter{ifeature}(iparam, :))
                
                xlabel(namesLocComb{combInd(icomb_real, 1)})
                ylabel(namesLocComb{combInd(icomb_real, 2)})
                end
            end % end of icomb
            
            set(findall(gcf, '-property', 'FontSize'), 'FontSize',18)
            sgtitle(sprintf('[%s] %s %s', namesType{itype}, namesFeature{ifeature}, namesParams_unit{iparam}), 'FontSize',20)
            
        end % end of iparam
    end % end of itype
end % end of ifeature

