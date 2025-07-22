% 
% 
% 
% if n == 1, iLocPlot = 1:nLoc; else, iLocPlot = 2:5;  end
% 
% % decide the axis range of colorbar
% data2D_toPlot = data2D(iLocPlot, :, :);
% caxisMin = min(data2D_toPlot(:));
% caxisMax = max(data2D_toPlot(:));
% 
% figure('Position', [0 0 1000 800])
% for iLoc = iLocPlot
%     subplot(3,3, subplot_locs(iLoc)), hold on
%     
%     % get data
%     data2D_ = flip(squeeze(data2D(iLoc , :, :))');
%     
%     %%%%%%%%%%%%%%%%
%     % work on this function !! %%%%
%     %%%%%%%%%%%%%%%%
%     
%     RCplot_2Dkernel(data2D_, xaxis_all{1}, xaxis_all{2}, ticks_ORI, ticks_SF, ticklabels_ORI, ticklabels_SF, caxisMin, caxisMax, namesLoc2D, plotPBorder, outline, xticks_border, yticks_border)
%     if iLoc == 2, xlabel('ORI (deg)'), ylabel('SF (cpd)'), end
%     
%     title(namesLoc2D{iLoc})
% end
% set(findall(gcf, '-property', 'FontSize'), 'FontSize',18)
