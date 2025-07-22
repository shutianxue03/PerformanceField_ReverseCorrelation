

namesTitles1 = {'Fovea', 'HM', 'Lower', 'Left', 'Fovea', 'Fovea'}; % before 'vs.', color-coded
namesTitles2 = {'Periphery', 'VM', 'Upper', 'Right', 'HM', 'VM'}; % after 'vs.', color-coded
namesTitles3 = {'', '', 'VM', 'HM', '', ''}; % in black

sz_text = 18;
sz_title = 20;
sz_ticks = 15;

yticks_ = {-.05:.05:.2, -.02:.03:.13};
ytext = {[.2, .18], [.13, .118]};
ylimit = {[-.06, .21], [-.026, .136]};

line_extra = [0, 1];

for itype = 1:ntypes
    figure('Position', [0 0 2000 800])
    
    for ifeature = 1:2
        xaxis = axis_tuning{ifeature};
        xaxis_itp = axis_itp_tuning{ifeature};
        nx = length(xaxis_itp);
        nameModelParams = namesParams_all{ifeature};
        nparams = length(namesParams_all);
        
        % ==== ylim is consistent within one trial type ====
        
        for icomb = 1:ncomb6
            if plotOne, subplot(2, ncomb6, plotInd(ifeature, icomb)),
            else, figure('Position', [0 200 600 500])
            end
            
            hold on, box on
            legends_all = {};
            
            color1 = colors_comb(combInd(icomb, 1),:);
            color2 = colors_comb(combInd(icomb, 2),:);
            
            for iline = 1:2 % 1=fovea/HM/LVM/HVM; 2=peri/VM/UVM/RHM
  
                % to locate the real location (given that only looking at dominant eyes)
                iLocComp = combInd(icomb, iline);

                % extract values and get ave & CI/SEM
                color = colors_comb(combInd(icomb, iline),:);
                
                % extract data
                marg_aveSubj = marg_aveSubj_all{ifeature, itype, iLocComp}.';
                marg_semSubj_neg = marg_semSubj_all{ifeature, itype, iLocComp}.';
                marg_semSubj_pos = marg_semSubj_neg;
                margPred_aveSubj = pred_aveSubj_all{ifeature, itype, icomb}(iline, :);
                margPred_semSubj_neg = pred_semSubj_all{ifeature, itype, icomb}(iline, :);
                margPred_semSubj_pos = margPred_semSubj_neg;
                margR2_aveSubj = R2_aveSubj_all{ifeature, itype, icomb}(iline);
                margR2_semSubj_neg = R2_semSubj_all{ifeature, itype, icomb}(iline);
                margR2_semSubj_pos = margR2_semSubj_neg;
                
                % extra lines
                yline(0, 'handlevisibility', 'off', 'linewidth', 1, 'color', [.7, .7, .7]);
                xline(line_extra(ifeature), 'handlevisibility', 'off', 'linewidth', 1, 'color', [.7, .7, .7]);
                
                % plot raw data
                errorbar(xaxis, marg_aveSubj, marg_semSubj_neg, marg_semSubj_pos, 'color', color, 'CapSize',0, 'linestyle', 'none')
                plot(xaxis, marg_aveSubj, 'o', 'linestyle', 'none', 'markersize', 1/sqrt(nx)*50, 'markerfacecolor', color, 'markeredgecolor','w')
                
                % plot fitting lines
                plot(xaxis_itp, margPred_aveSubj, '-', 'color', color, 'linewidth', 2)
                patch([xaxis_itp, flip(xaxis_itp)], [margPred_aveSubj-margPred_semSubj_neg, flip(margPred_aveSubj+margPred_semSubj_neg)].',  color, 'FaceAlpha', .3, 'linestyle', 'none')
                
                %  print R2
                if nsubj == 1, text(axis_tuning{ifeature}(1), ytext{ifeature}(iline), sprintf('R^2 = %.3f', margR2_aveSubj), 'color', color, 'fontsize', sz_text)
                else, text(axis_tuning{ifeature}(1), ytext{ifeature}(iline), sprintf('R^2 = %.3f (+-%.2f)', margR2_aveSubj, margR2_semSubj_neg), 'color', color, 'fontsize', sz_text)
                end
                % plot tuning curve of gabor and calculate the corr coeff
                %                 marg_gabor = reshape(marg_gabor, length(marg_gabor), 1);
                % %                 plot(xaxis, marg_gabor, 'g-')
                %                 [r, p] = corr([marg_gabor, margPred_ave]);
                %                 string_stars = getString_starts( p(1,2));
                %                 text_corr = sprintf('r=%.2f%s',r(1,2), string_stars);
                %                 text(xticks_{ifeature}(1), ymin+.015*iline, text_corr, 'color', color, 'fontsize', sz_text)
                %
            end
            
            % ==================
            % combare two lines (nB>1)
            %             if nB>1
            %                 y = squeeze(marg{ifeature}(:, itype, icomb, :, :));
            %                 for ix = 1:nx
            %                     if icomb == 4 % test if LEFT is different from RIGHT
            %                         [h,p] = ttest(squeeze(y(:, 1, ix)), squeeze(y(:, 2, ix)));
            %                     else % for other comparisons, also test direction
            %                         [h,p] = ttest(squeeze(y(:, 1, ix)), squeeze(y(:, 2, ix)), 'Tail', 'right');
            %                     end
            %                     p = p*nx;
            %                     if p<.001, plot(xaxis_itp(ix), ylimit{ifeature}(2)+.01, 'k*'), end
            %                 end
            %             end
            % ==================
            
            % xaxis
            xlim(axisLim{ifeature})
            xlabel(namesFeature{ifeature}, 'fontsize', sz_ticks)
            xticks(axisTicks_tuning{ifeature})
            xticklabels(axisTL_tuning{ifeature})
            
            % yaxis
            ylim(ylimit{ifeature})
            yticks(yticks_{ifeature})
%             if icomb == 1, ylabel([titleF{ifeature}, ' (a.u.)'], 'fontsize', sz_ticks), end
            
            ax = gca;
            ax.FontSize = sz_ticks;
            
            % make colorful titles
            a = ['\color[rgb]{', sprintf('%1.2f, %1.2f, %1.2f', color1), '} ',  namesTitles1{icomb}];
            b = '\color[rgb]{0,0,0} vs. ';
            c = ['\color[rgb]{', sprintf('%1.2f, %1.2f, %1.2f', color2), '}',  namesTitles2{icomb}];
            d = [' \color[rgb]{', sprintf('%1.2f, %1.2f, %1.2f', [0 0 0]), '}',  namesTitles3{icomb}];
            
            set(findall(gcf, '-property', 'linewidth'), 'linewidth',1.5)
            
            if plotOne
                if ifeature==1, title([a,b,c,d], 'fontsize', sz_title, 'interpreter', 'tex'), end
            else
                title([a,b,c,d], 'fontsize', sz_title, 'interpreter', 'tex')
                saveas(gcf, sprintf('publishedPDFs/fig/SP_%s_comb%d.jpg', kernelNames{ifeature}, icomb))
            end
        end % end of icomb
    end % end of ifeature
    
%     make_sgtitle(sprintf('[%s-%s] Marginalized values ', titleF{ifeature}, namesType{itype}), subjName, nB, nsubj, nAllTrials)
end % end of itype





