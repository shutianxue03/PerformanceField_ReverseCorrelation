
subplot_locs = [5,4,2,6,8];
%% extract titration data
fileNameT = sprintf('%s%s_stair', folderName, subjName);
fileNameT_dir = dir([fileNameT,'*']); % find out all possible files
if isempty(fileNameT_dir), error('ALERT: titration file of this subject is absent !!')
else, load(fileNameT_dir.name), stairsT = record.stairs_all;
end
clear record

%% extract Qstair data
fileNameQ = sprintf('%s%s_qstair', folderName, subjName);
fileNameQ_dir = dir([fileNameQ,'*']); % find out all possible files
load(fileNameQ_dir.name)
stairsQ = record.stairs_all;
clear record

%%
figure('Position',[0 200 800 800])
for iLoc = 1:design.nPrsLoc
    subplot(3,3,subplot_locs(iLoc)), hold on
    
    x = stairsQ{iLoc}.x(2:end);
    
    % plot the two endpoints of the titration
    plot(length(x)+5,  stairsT{iLoc,1}.x(end), 'ro') 
    plot(length(x)+5,  stairsT{iLoc,2}.x(end), 'bo')
    
    % plot the staircase of the qStair
    plot(x, 'k', 'handlevisibility', 'off')
    
    % plot the endpoint of the qStair
    plot(length(x), x(end), '*')
    
    if iLoc == 2, xlabel('trial #'), ylabel('target cst (%)'), legend({'titration1','titration2','qstair', 'last thresh'}, 'Location', 'best'),end
    
    title(sprintf('thresh = %.2f', x(end)))
end
set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)
sgtitle(sprintf('%s - quick staircase (%d trials)', subjName, record.params.design.nTrialsPerBlock), 'fontsize', 20);

%% obtain and plot PDF
plotFlag = 1;

figure('Position',[800 200 800 800])
flattenFactor = 1;

stairs_inUse = stairsQ;

for iLoc = 1:design.nPrsLoc
    ymax = .015;
    x = stairs_inUse{iLoc}.priorAlphaRange;
    getPDF
    plot(xlog, stairsT{iLoc,1}.pdf, 'r--')
    plot(xlog, stairsT{iLoc,2}.pdf, 'b--')
    if iLoc == 2, legend('qstair pdf', 'prior1', 'prior2'), end
end


set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)
sgtitle(sprintf('%s - quick staircase (%d trials)', subjName, recordQ.params.design.nTrialsPerBlock), 'fontsize', 20);

%     compareThresh