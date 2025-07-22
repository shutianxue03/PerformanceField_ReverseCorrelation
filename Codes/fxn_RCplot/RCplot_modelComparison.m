
% plot the group averaged AICc


for ifeature = 2
    example = delta_allSubj{1};
    nmodels = length(example{1,ifeature,1});

    for icomb = 1:4
        figure('Position', [0 200 40*nmodels 1200])
        for itype = 1:ntypes
            
            delta_allSubj2 = nan(nsubj, nmodels);
            for isubj = 1:nsubj
                aa = delta_allSubj{isubj};
                delta_allSubj2(isubj, :) = aa{icomb, ifeature, itype};
            end
            
            %% plot
            % decide the format of the x tick labels (0=short, M%d, 1=long, plus parameter comb)
%             if itype == ntypes, flagLongXTicks = 1; else, flagLongXTicks = 0; end
            flagLongXTicks = 1;
            subplot(ntypes, 1, itype), hold on
            
            % =========
            RCplot_AICc
            % =========
            title(sprintf('[%s-%s] %s vs. %s', namesFeature{ifeature}, namesType{itype}, namesLocComb{combInd(icomb, [1,2])}))
            set(findall(gcf, '-property', 'FontSize'), 'FontSize',12)
            set(findall(gcf, '-property', 'linewidth'), 'linewidth', .5)
            
        end
    end
end

