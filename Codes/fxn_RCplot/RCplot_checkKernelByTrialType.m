
% load(energyMatDir.name)
clear patch

figure('Position', [0 200 1200 300])
lineCol = {'r', 'b', 'k'};
faceAlpha = 0.2; % the transparency of the shade

fileDir_boot = dir(fileName_boot);

if nsubj > 1
    kernel_allLoc = squeeze(mean(kernel_allSubj, 1));
    kernel_SEM = squeeze(std(kernel_allSubj, [], 1))/sqrt(nsubj);
end

for iLoc = 1:nLoc
    
    subplot(1, nLoc, iLoc), hold on
    
    for itype = 1:ntypes   % prs, abs, both
        kernel_ = squeeze(kernel_allLoc(itype, iLoc, :));
        if nsubj>1, kernelSEM_ = squeeze(kernel_SEM(itype, iLoc, :));end
        
        % plot the kernel
        plot(filterSF_all_log, kernel_,[lineCol{itype}, '-'])
        
        % plot SEM
        if nsubj == 1 
            if ~isempty(fileDir_boot), load(fileDir_boot.name)
            ub = squeeze(kernel_68CI(1,itype,iLoc,:));
            lb = squeeze(kernel_68CI(2,itype,iLoc,:));
            end
        else
            ub = kernel_ + kernelSEM_;
            lb = kernel_ - kernelSEM_;
        end
        x = [filterSF_all_log, flip(filterSF_all_log)];
        y = [lb; flip(ub)];
        patch('Faces',1:2*nfilters,'Vertices', [x', y], 'FaceColor', lineCol{itype}, 'FaceAlpha', faceAlpha,'EdgeColor','none', 'handlevisibility', 'off');
        
%         plot idvd data (given multiple subj)
%         if nsubj > 1, plot(repmat(filterSF_all_log, nsubj, 1), squeeze(kernel_allSubj(:,itype, iLoc, :)), '.', 'color', lineCol{itype}), end
        
    end
    
    plot(log2([.8, 4.8]), [0,0], 'color', [.5, .5, .5])
    plot(log2([2,2]), [-.2, 1], 'color',[.5,.5,.5])
    
    xticks(log2([1,2,4])), xticklabels([1,2,4])
    xlim(log2([.8,4.8]))
    title(locNames{iLoc})
    ylim(ylimit)
    if iLoc == 1, xlabel('SF channel (cpd)'), ylabel('kernel (a.u.)'), legend('tgt-prs trials', 'tgt-abs trials', 'all trials', 'location', 'best'), end
end

set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)
if nsubj == 1, sgtitle(sprintf('%s: %d trials per loc', subjName, nAllTrials), 'FontSize',25)
else, sgtitle(sprintf('n = %d (%d +- %d)', nsubj, round(mean(nAllTrials)), round(std(nAllTrials)/sqrt(nsubj))), 'FontSize',25)
end

