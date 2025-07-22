
iplotInd = [1,2,3;4,5,6];

for itype = 1:ntypes
    figure('Position', [0 0 1200 400])
    for icomb = 6:8
        %         figure('Position', [0 550 800 300])
        %         kernel2D1 = squeeze(kernels2D_perComb(itype, icomb, :, :));
        %         kernel2D2 = squeeze(kernels2D_perComb2(itype, icomb, :, :));
        %         subplot(1,3,1), imagesc(kernel2D1), title('By concat'), axis square
        %         subplot(1,3,2), imagesc(kernel2D2), title('By ave'), axis square
        %         subplot(1,3,3), imagesc(kernel2D1-kernel2D2); colorbar,  title('By concat minus By ave'), axis square
        %         sgtitle(sprintf('[%s] %s', namesType{itype}, namesLocComb{icomb}))
        %         set(findall(gcf, '-property', 'FontSize'), 'FontSize',18)
        
        %%%%%%
        %    ORI   %
        %%%%%%
        margOri1 = squeeze(margOri_perComb(itype, :, icomb));
        margOri2 = margOri_perComb2{icomb}(itype, :);
        
        subplot(2,3,iplotInd(1, icomb-5)), hold on
        plot(axis_tuning{1}, margOri1, ['-', colorsType{itype}])
        plot(axis_tuning{1}, margOri2, ['--', colorsType{itype}])
        xline(0);
        xticks(-90:45:90), xticklabels(-90:45:90)
        xlim([-90, 90])
        if icomb == 6, ylabel('ORI kernels'), end
        title(namesLocComb{icomb})
        %         subplot(2,2,3), hold on, axis square, box on
        %         plot(margOri1, margOri2, 'o')
        %         plot([-.05, .1], [-.05, .1], 'k-'), xlim([-.05,.1]), ylim([-.05,.1])
        %         xlabel('By concat')
        %         ylabel('By averaging')
        
        %%%%%%
        %    SF   %
        %%%%%%
        margSF1 = squeeze(margSF_perComb(itype, :, icomb));
        margSF2 = margSF_perComb2{icomb}(itype, :);
        
        subplot(2,3,iplotInd(2, icomb-5)), hold on
        plot(axis_tuning{2}, margSF1, ['-', colorsType{itype}])
        plot(axis_tuning{2}, margSF2, ['--', colorsType{itype}])
        xticks([0,1,2]), xticklabels([1,2,4])
        xlim([0,2])
        xline(1);
        if icomb == 6
            ylabel('SF kernels')
            legend({'by concate', 'by averaging'}, 'Location', 'best')
        end
        
        %         subplot(2,2,4), hold on, axis square, box on
        %         plot(margSF1, margSF2, 'o')
        %         plot([-.05, .1], [-.05, .1], 'k-'), xlim([-.05,.1]), ylim([-.05,.1])
        %         sgtitle(sprintf('[%s] %s', namesType{itype}, namesLocComb{icomb}))
        
    end % icomb
    set(findall(gcf, '-property', 'FontSize'), 'FontSize',16)
    set(findall(gcf, '-property', 'linewidth'), 'linewidth',1.5)
    sgtitle(namesType{itype})
end