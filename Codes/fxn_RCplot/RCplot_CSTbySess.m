

hold on, grid on

for iLoc = 1:nLoc
    if iLoc== 4, plot(log10(cst_perSess_perLoc(:,iLoc)), '--', 'color',colors_comb(iLoc, :)),
    else, plot(log10(cst_perSess_perLoc(:,iLoc)), 'color',colors_comb(iLoc, :)), end
end

yticks_log = quantile(log10(cst_perSess_perLoc(:)), 0:.5:1);
yticks(yticks_log)
yticklabels(round(10.^yticks_log, 2))
xlabel('Session #'), ylabel('Gabor CST')
if nSess > 1, xlim([.5,nSess + .5])
    if nSess>2, xticks([1, round(nSess/2),nSess]), xticklabels([1, round(nSess/2),nSess]),else, xticks([1,2]), end
else, xlim([.5,2]), xticks(1)
end
%     xlim([.5, nSess+.5])
legend(namesLoc2D, 'location', 'eastoutside')