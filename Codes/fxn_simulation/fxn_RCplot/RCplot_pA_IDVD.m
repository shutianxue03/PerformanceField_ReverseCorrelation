
% pA_prs_perLoc/pA_abs_perLoc are calculated in 'RC_pA.m'
% 1st row: resp=0 given PRS/ABS
% 2nd row: inconsistent resp given PRS/ABS
% 3rd row: resp=1 given PRS/ABS



%% plot pA across loc, considering signal existence
%
figure, hold on
buff = .1;
pAprsColors = [1,2/3,2/3; 1,1/3,1/3; 1,0,0;  ];
pAabsColors = [0,0,1;  1/3,1/3,1; 2/3,2/3,1];

% gabor-PRS
for n = [3,2,1]
    pA_perSess_perLoc_ave = squeeze(mean(pA_perSess_perLoc(n, :, :)));
    pA_perSess_perLoc_sd = squeeze(std(pA_perSess_perLoc(n, :, :)))/sqrt(nSess);
    
    if n==2, lineStyle = '--'; else, lineStyle = '-'; end
    plot(1:nLoc8, pA_prs_perLoc(:, n), 'o', 'Linestyle', lineStyle, 'color', pAprsColors(n, :))
end

% gabor-ABS
for n = 1:3
    if n==2, lineStyle = '--'; else, lineStyle = '-'; end
    plot(1:nLoc8, pA_abs_perLoc(n, :), 'o', 'Linestyle', lineStyle, 'color', [n/3-buff, n/3-buff, 1])
end

legend({'both resp=prs', 'inconsistent resp', 'both resp=abs', 'both resp=abs', 'inconsistent resp', 'both resp=prs'}, 'location', 'eastoutside')
xticks(1:5), xticklabels(namesLoc2D), xlim([.5, 5.5])
ylabel('Proportion')
title('Proportion of Gabor-\color[rgb]{1,0,0}PRS \color[rgb]{0,0,0}/\color[rgb]{0,0,1}ABS\color[rgb]{0,0,0} trials')

set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)
set(findall(gcf, '-property', 'linewidth'), 'linewidth',1.5)


