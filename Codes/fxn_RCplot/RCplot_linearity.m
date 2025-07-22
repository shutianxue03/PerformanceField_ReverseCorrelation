
%%
ylabels_pC = {'rHit', 'rCR', 'pC'};
ncomp = 3; % 1: x=pA, y=pC, 2: x=cst, y=pC; 3: x=cst, y=pA
colorUniks = rand(30, 3); % one color for each cst level

%%
for itype = 1:ntypes
    for iLoc = 1:nLoc5
        [cst_unik_all, ~, iC_cst] = unique(cst_perSess_perLoc(:, iLoc)); % unique contrasts used at this loc
        ncst_unik_all = length(cst_unik_all); % number of unique contrasts used at this loc
        
        nblocks_perCST = nan(1, ncst_unik_all); % number of blocks per contrast
        pA_perCST = nblocks_perCST;
        pC_perCST = nblocks_perCST;
        
        for icst_unik = 1:ncst_unik_all
            cst_unik = cst_unik_all(icst_unik);
            nblocks_perCST(icst_unik) = sum(cst_unik==cst_perSess_perLoc(:, iLoc)); % number of blocks per contrast
            pA_perCST(icst_unik) = mean(pA_perSess_perLoc(itype, cst_unik==cst_perSess_perLoc(:, iLoc), iLoc));
            pC_perCST(icst_unik) = mean(pC_perSess_perLoc(itype, cst_unik==cst_perSess_perLoc(:, iLoc), iLoc));
        end
        
        nblocks_unik_perLoc{iLoc} = nblocks_perCST; % unik stands for unique
        cst_unik_perLoc{iLoc} = cst_unik_all;
        pA_unik_perLoc{iLoc} = pA_perCST;
        pC_unik_perLoc{iLoc} = pC_perCST;
        
    end
    
    figure('Position', [0 0 1500 ncomp*300])
    for icomp = 1:ncomp
        switch icomp
            case 1 % pA vs. pC
                x = squeeze(pA_perSess_perLoc(itype, :, :));
                y = squeeze(pC_perSess_perLoc(itype, :, :));
                x_unik = pA_unik_perLoc;
                y_unik = pC_unik_perLoc;
                xlabel_ = 'pA';
                ylabel_ = ylabels_pC{itype};
            case 2 % CST vs. pC
                x = log10(cst_perSess_perLoc);
                y = squeeze(pC_perSess_perLoc(itype, :, :));
                x_unik = cst_unik_perLoc;
                y_unik = pC_unik_perLoc;
                xlabel_ = 'CST';
                ylabel_ = ylabels_pC{itype};
            case 3 % CST vs. pA
                x = log10(cst_perSess_perLoc);
                y = squeeze(pA_perSess_perLoc(itype, :, :));
                x_unik = cst_unik_perLoc;
                y_unik = pA_unik_perLoc;
                xlabel_ = 'CST';
                ylabel_ = 'pA';
        end
        
        % get xlimit and ylimit
        xmin = min(x(:));
        xmax = max(x(:));
        ymin = min(y(:));
        ymax = max(y(:));
        
        for iLoc = 1:nLoc5
            subplot(ncomp, nLoc5, (icomp-1)*nLoc5+iLoc), hold on, grid on, box on
            cst_unik_all = cst_unik_perLoc{iLoc};
            nblocks_perCST = nblocks_unik_perLoc{iLoc};
            ncst_unik_all = length(nblocks_perCST);
            
            for icst_unik = 1:ncst_unik_all
                % idvd block (NOT averaged across cst)
                ind = cst_unik_all(icst_unik) == cst_perSess_perLoc(:, iLoc);
                if size(colorUniks,1) < icst_unik, error('ALERT: not enough color options!'), end
                
                % plot data points at EACH CST level
                if plotFlag_unik
                    plot(x(ind, iLoc), y(ind, iLoc),  'x', 'color', colorUniks(icst_unik, :))
                end
                % plot data points AVERAGED across CST levels
                markerSize = nblocks_perCST(icst_unik)/sum(nblocks_perCST)*10+10;
                if icomp ~= 1
                    x_unik_perLoc = log10(x_unik{iLoc});
                else
                    x_unik_perLoc = x_unik{iLoc};
                end
                
                if plotFlag_unik
                    plot(x_unik_perLoc(icst_unik), y_unik{iLoc}(icst_unik), 'o', 'color', colorUniks(icst_unik, :), 'markersize', markerSize)
                else
                    plot(x_unik_perLoc(icst_unik), y_unik{iLoc}(icst_unik), 'o', 'color', colors_comb(iLoc, :), 'markersize', markerSize)
                end
            end
            
            % linear regression
            if plotFlag_unik
                x_ = x(:, iLoc);
                y_ = y(:, iLoc);
            else
                x_ = reshape(x_unik_perLoc, length(x_unik_perLoc), 1);
                y_ = reshape(y_unik{iLoc}, length(x_unik_perLoc), 1);
            end
            
            lm = polyfit(x_, y_, 1);
            [r,p] = corr(x_, y_);
            yfit = polyval(lm, x_);
            plot(x_, yfit, 'color', colors_comb(iLoc, :), 'handlevisibility', 'off')
            R2 = getR2(y_, yfit);
            ss = getString_starts(p);
            
            if icomp == 1, title(sprintf('%s\nr = %.2f%s', namesLoc2D{iLoc}, r, ss))
            else, title(sprintf('r = %.2f%s', r, ss))
            end
            
            if (icomp == 1) || (iLoc == 1)
                xlabel(xlabel_)
                ylabel(ylabel_)
            end
            xlim([xmin, xmax])
            ylim([ymin, ymax])
            
            if icomp ~= 1, xticks(xmin:.1:xmax), xticklabels(round(10.^(xmin:.1:xmax)*100)), end
            %         if icomp == 1, title(namesLoc2D{iLoc}), end
        end % end of iLoc
        
    end % end of icomp
    
    set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)
    set(findall(gcf, '-property', 'linewidth'), 'linewidth',1)
    
    make_sgtitle(sprintf('[%s] CST vs. pA vs. pC', namesType{itype}), subjName, nB, nsubj, nAllTrials)
end % end of itype


%% Linear mixed model
LocInd = repmat(1:5, ntypes, 1, nblocks/nLoc5);
LocInd = reshape(LocInd, [ntypes, nSess, nLoc5]);
typeInd = repmat(1:3, nblocks/nLoc5, 1, nLoc5);
typeInd = reshape(typeInd, [ntypes, nSess, nLoc5]);

% predict pA by itype and iLoc
printLMM({pA_perSess_perLoc, typeInd, LocInd},  {'pA','Type','Loc'}, 'pA~Type+Loc')

% predict pC by itype and iLoc
printLMM({pC_perSess_perLoc, typeInd, LocInd},  {'pC','Type','Loc'}, 'pC~Type+Loc')

