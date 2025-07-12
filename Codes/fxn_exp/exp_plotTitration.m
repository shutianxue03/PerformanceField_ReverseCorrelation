
folderName = ['Data/', subjName, '/'];
ratioSelected = .7; 
isubj = 1; 
thresh_SEM = zeros(1,5); 
plotFlag = 1;
draw_mode = 1; 
plotErrorBar = 0;

%% get PDF and save it separatly
figure('Position',[0 200 800 800])
x = record.params.stairParams.alphaRange;
nx = length(x);
lastPosterior_allLoc = nan(5, nx);
mean_allLoc = nan(5, 1);
std_allLoc = nan(5, 1);
flattenFactor = 1;
for iLoc = params.design.loc_list %CHANGED: 1:5locs
    ymax = .01;
    %------%
    getPDF
    %------%
    lastPosterior_allLoc(iLoc,:) = pdf_ave_flatten;
    mean_allLoc(iLoc) = mean_comb;
    std_allLoc(iLoc) = sd_flat;
end

set(findall(gcf, '-property', 'FontSize'), 'FontSize',12)
% sgtitle(sprintf('S%d %s (%d sessions)', isubj, subjName, nSess), 'fontsize', 20);
save([folderName,subjName, '_PDF'], 'x', 'lastPosterior_allLoc', 'mean_allLoc', 'std_allLoc')
saveas(gcf, sprintf('%s/FIG_ttt_PDF.jpg', folderName))

%% plot staircase and SDT data
drawEachSubj

set(findall(gcf, '-property', 'FontSize'), 'FontSize',10)

% subplot(3,3,9)
% plot(x, lastPosterior_allLoc)
% xlim([0,0.3])
% legend({'center', 'left', 'UVM', 'right', 'LVM'}, 'Location', 'best')
% title('all flattened PDF')
% set(findall(gcf, '-property', 'LineWidth'), 'LineWidth',2)
% set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)
