

for iModelB = 1:nModelsB
    for iModelA = 1:nModelsA
        figure('Position', [3e3 0 1000 500])
        
        maxChannel_allSubj_ = squeeze(NOM_maxChannel_allSubj(:, iiLoc, iModelA, iModelB, :, :));
        [maxChannel_ave, ~, ~, maxChannel_sem] = getCI(maxChannel_allSubj_, 2, 1);
        maxChannel_sd = maxChannel_sem * sqrt(nsubj);
        
        for iFeature = 1:2
            %%%%%%%%%
            % bar
            %%%%%%%%%
            subplot(3,2, iFeature), hold on,
            for iPRS = 1:2
                errorbar(iPRS, maxChannel_ave(iPRS, iFeature), maxChannel_sd(iPRS, iFeature), 'ko', 'CapSize', 0)
                xlim([0, 3])
                xticks(1:2)
                xticklabels({'PRS', 'ABS'})
                
                if iFeature == 1, yline(0, 'r-'); ylim([-40, 40]), title('ORI')
                else, yline(1, 'r-'); ylim([.5, 1.5]), yticks(.5:.25:1.5), yticklabels(round(2.^(.5:.25:1.5), 2)), title('SF')
                end
            end % iPRS
            
            %%%%%%%%%
            % idvd data 
            %%%%%%%%%
            for iPRS = 1:2
                subplot(3,2, iFeature+iPRS*2), hold on,
                plot(1:nsubj, maxChannel_allSubj_(:, iPRS, iFeature), 'o')
                if iFeature == 1, yline(0, 'r-'); ylim([-40, 40])
                else, yline(1, 'r-'); ylim([.5, 1.5]), yticks(.5:.25:1.5), yticklabels(round(2.^(.5:.25:1.5), 2))
                end
                title(namesType{iPRS})
                xticks(1:nsubj)
                xticklabels(subjList)
                xlim([0, nsubj+1])
            end
            
        end % iFeature
        sgtitle(sprintf('[%s] A%dB%d\nError bar is SD', namesLocComb{iLoc}, iModelA, iModelB))
    end % iModelA
    
    %     if iModelA==1, saveas(gcf, sprintf('talk2022/fig/model/imax.jpg')), end
end % iModelB

