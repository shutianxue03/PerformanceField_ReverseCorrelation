

function RCplot_2D(data2D, LocFlag, nB, title_)
SX_RC1_setting
pCriterion = .001;
data2D_ave = getCI(data2D, 1, 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5 loc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for itype = 1:ntypes
%     data2D_ave_perType = squeeze(data2D_ave(itype, :, :, :));
%     data2D_perType = squeeze(data2D(:, itype, :, :, :));
%
%     for n=1:LocFlag
%         if n == 1, iLocPlot = 1:nLoc5; else, iLocPlot = 2:5;  end
%         % decide the axis range of colorbar
%         data2D_toPlot_ave = data2D_ave_perType(iLocPlot, :, :);
%
%         caxisLim = [min(data2D_toPlot_ave(:)), max(data2D_toPlot_ave(:))];
%         figure('Position', [0 0 1000 800])
%         for iLoc = iLocPlot
%             subplot(3,3, subplot_locs(iLoc)), hold on
%             % get data
%             data2D_ave_perLoc = squeeze(data2D_ave_perType(iLoc , :, :));
%             data2D_perLoc = squeeze(data2D_perType(:, iLoc, :, :));
%
%             % get outline
%             if nB > 1
%                 outline = getOutline(nORI, nSF, data2D_perLoc, nB, pCriterion);
%                 if strcomp(title_, 'p values'), outline = data2D_ave_ < pCriterion; end
%             else, outline = [];
%             end
%
%             % plot
%             RCplot_2Dkernel(data2D_ave_perLoc', outline', caxisLim)
%
%             if iLoc == 2, xlabel('ORI (deg)'), ylabel('SF (cpd)'), end
%             title(namesLoc2D{iLoc})
%         end % end of iLoc
%
%         set(findall(gcf, '-property', 'FontSize'), 'FontSize',18)
%         make_sgtitle(sprintf('[%s] %s', namesType{itype}, title_) , subjName, nB, nsubj, nAllTrials_allSubj)
%
%     end % end of n = 1:LocFlag
% end % end of itype
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pre3 and comparison of combined loc %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for itype = 1:ntypes
    figure('Position', [0 0 ncomb6*300 1000])
%     for icomb = 1:ncomb6
        data2D1_ave = squeeze(mean(data2D_ave(:, 1, itype, :, :), 1));
        data2D2_ave = squeeze(mean(data2D_ave(:, 2, itype, :, :), 1));
        
        data2D1 = squeeze(data2D_ave(:, 1, itype, :, :));
        data2D2 = squeeze(data2D_ave(:, 2, itype, :, :));
%         data2D2 = squeeze(data2D(:, itype, combInd(icomb, 2), :, :));
        
        %%%%%%%%%%%%%%%%%%%%%%
        %   1. 2D kernels of Fov/HM/LVM/Left   %
        %%%%%%%%%%%%%%%%%%%%%%
%         if nB>1, 
            outline = getOutline(nORI, nSF, data2D1 - data2D2, nB, pCriterion);
%         else, outline = [];
%         end
        caxisLim = [min(data2D1_ave(:)), max(data2D1_ave(:))];
        if abs(caxisLim(1) - caxisLim(2)) <= eps, caxisLim = caxisLim +[-.001, .001]; end
        
        subplot(3,1,1), hold on
        RCplot_2Dkernel(data2D1_ave', outline', caxisLim)
        title(namesLocComb{combInd(icomb, 1)})
        
        %%%%%%%%%%%%%%%%%%%%%%%
        %   2. 2D kernels of Peri/VM/UVM/Right   %
        %%%%%%%%%%%%%%%%%%%%%%%
        caxisLim = [min(data2D2_ave(:)), max(data2D2_ave(:))];
        if abs(caxisLim(1) - caxisLim(2)) <= eps, caxisLim = [caxisLim(1) -.001,  caxisLim(2) +.001];end
        if nB>1, outline = getOutline(nORI, nSF, data2D2, nB, pCriterion);
        else, outline = [];
        end
        subplot(3,ncomb6, icomb+ncomb6), hold on
        RCplot_2Dkernel(data2D2_ave', outline', caxisLim)
        title(namesLocComb{combInd(icomb, 2)})
        
        %%%%%%%%%%%
        %       3. diff            %
        %%%%%%%%%%%
        data2D1_allB = squeeze(data2D(:, itype, combInd(icomb, 1), :, :));
        data2D2_allB = squeeze(data2D(:, itype, combInd(icomb, 2), :, :));
        data2D_diff_allB = data2D1_allB - data2D2_allB;
        
        if nB>1, data2D_diff_ave = squeeze(mean(data2D_diff_allB, 1));
        else
            if nsubj==1, data2D_diff_ave = data2D_diff_allB;
            else, data2D_diff_ave = squeeze(mean(data2D_diff_allB, 1));
            end
        end
        caxisLim = [min(data2D_diff_ave(:)), max(data2D_diff_ave(:))];
        
        subplot(3,ncomb6, icomb+ncomb6*2), hold on
        if nB>1
            % get outline
            outline = getOutline(nORI, nSF, data2D_diff_allB, nB, pCriterion);
            RCplot_2Dkernel(data2D_diff_ave', outline', caxisLim)
        else
            RCplot_2Dkernel(data2D_diff_ave', [], caxisLim)
        end
        
        title(sprintf('%s minus %s', namesLocComb{combInd(icomb, 1)}, namesLocComb{combInd(icomb, 2)}))
        
%     end % end of icomb
    set(findall(gcf, '-property', 'FontSize'), 'FontSize',18)
    make_sgtitle(sprintf('[%s] %s', namesType{itype}, title_) , subjName, nB, nsubj, nAllTrials)
    
end % end of itypes



