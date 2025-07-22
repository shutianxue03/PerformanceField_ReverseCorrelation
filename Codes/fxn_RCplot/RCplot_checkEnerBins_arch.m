
%% this is arch because only SF channels were considered

slopes = nan(nfilters, nLoc, 2);
nfilter_s = length(filter_s);


%% plot hit/FA rate as a fxn of energy (unnormalized)
figure('Position', [0 300 nfilter_s*400 600])

if nsubj>1
    x_prs_ = squeeze(mean(x_prs_allSubj, 1)); % x_prs_allSubj: nsubj x nLoc x nfilters x nbins_e
    x_prs_SEM = squeeze(std(x_prs_allSubj, [], 1));
    x_abs_ = squeeze(mean(x_abs_allSubj, 1));
    x_abs_SEM = squeeze(std(x_abs_allSubj, [], 1));
    rHit_ = squeeze(mean(rHit_allSubj, 1));
    rFA_ = squeeze(mean(rFA_allSubj, 1));
    rHit_SEM = squeeze(std(rHit_allSubj, [], 1));
    rFA_SEM = squeeze(std(rFA_allSubj, [], 1));
else
    x_prs_ = x_prs;
    x_prs_SEM = [];
    x_abs_ = x_abs;
    x_abs_SEM = [];
    rHit_ = rHit;
    rFA_ = rFA;
    rHit_SEM = [];
    rFA_SEM = [];
end

for ifilter = 1:nfilters
    for n=1:2 % prs and abs % x_prs and x_abs are from function checkEnergyByTrialType
        if n==1, x = squeeze(x_prs_(:, ifilter, :)).'; y = squeeze(rHit(:, ifilter, :)).'; xlimit = [0,5]; ylabel_ = 'hit rate';
            if nsubj>1
                x_anova = x_prs_allSubj; 
                x_SEM = squeeze(x_prs_SEM(:, ifilter, :));
                y_SEM = squeeze(rHit_SEM(:, ifilter, :));
                text_y = .25;
            end
            
        else,     x = squeeze(x_abs_(:, ifilter, :)).'; y = squeeze(rFA(:, ifilter, :)).'; xlimit = [0,2]; ylabel_ = 'FA rate';
            if nsubj>1
                x_anova = x_abs_allSubj; 
                x_SEM = squeeze(x_abs_SEM(:, ifilter, :)); 
                y_SEM = squeeze(rFA_SEM(:, ifilter, :));
                text_y = .75;
            end
        end
        
        % get slope
        for iLoc = 1:nLoc
            beta = polyfit(x(:, iLoc),y(:, iLoc),1);
            slopes(ifilter, iLoc, n) = beta(1);
        end
        
        % plot the hit/FA rate vs. binned energy
        xmin = floor(min(x(:)));
        xmax = ceil(max(x(:)));
        iplot = find(ifilter == filter_s);
        if iplot
            % ANOVA 2
            data_anova = [];
            ind_subj = [];
            ind_e = [];
            ind_l = [];
            
            for isubj = 1:nsubj
                ind_subj = [ind_subj; ones(nbins_e-1, nLoc) * isubj];
                ind_e = [ind_e; repmat((1:(nbins_e-1))', 1, nLoc)];
                ind_l = [ind_l; repmat(1:nLoc, nbins_e-1, 1)];
                data_anova = [data_anova; squeeze(x_anova(isubj,:,ifilter,:))'];
            end
            
            varNames = {'Loc', 'energy', 'Interaction'};
            [~, tbl, ~] = anova2(data_anova, nsubj, 'off');
            df2 = tbl{5,3};
            anova_text = [];
            for aa = 1:length(varNames)
                df1 = tbl{aa+1, 3};
                F = tbl{aa+1, 5};
                p =  tbl{aa+1, 6};
                
                anova_text = [anova_text, sprintf('%s: F(%d,%d) = %.3f, p = %.3f\n', varNames{aa}, df1, df2, F,p)];
            end
            
            % plot
            subplot(2, nfilter_s, iplot+(n-1)*nfilter_s), hold on, grid on
            for iLoc = 1:nLoc
                plot(x(:,iLoc), y(:,iLoc), '.-', 'color', colors_loc5(iLoc, :))
                if nsubj>1
%                     errorbar(x(:, iLoc), y(:, iLoc), x_SEM(iLoc, :), 'horizontal', 'color', colors_loc5(iLoc, :), 'handlevisibility', 'off')
                    errorbar(x(:, iLoc), y(:, iLoc), y_SEM(iLoc, :), 'vertical', 'color', colors_loc5(iLoc, :), 'handlevisibility', 'off')
                end
            end
            xlim([xmin, xmax]), ylim([0 1]), yticks(0:.25:1), title(sprintf('filter = %.2f cpd\n', filterSF_all(ifilter)))
            if iplot  == 1, ylabel(ylabel_), if n==1, xlabel('energy'), legend(locNames, 'location', 'southeast'), end, end
            text(xmin, text_y, anova_text)
        end
    end
end

set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)
sgtitle(sprintf('Behavior as a function of energy (%d bins)', nbins_e-1), 'FontSize', 20)

%% plot linear slopes
if nsubj==1
    figure('position', [0 0 1500 300])
    for iLoc = 1:nLoc
        subplot(1, nLoc, iLoc), hold on
        plot(filterSF_all_log, squeeze(slopes(:, iLoc, 1)).', 'r')
        plot(filterSF_all_log, squeeze(slopes(:, iLoc, 2)).', 'b')
        plot(filterSF_all_log(filter_s), squeeze(slopes(filter_s, iLoc, 1)).', 'ro')
        plot(filterSF_all_log(filter_s), squeeze(slopes(filter_s, iLoc, 2)).', 'bo')
        xline(log2(2), 'k', 'color',[.5,.5,.5]);
        yline(0, 'k', 'color',[.5,.5,.5]);
        
        if iLoc == 1, legend({'tgt-prs', 'tgt-abs'}), xlabel('SF channels'), ylabel('slope'), end
        xticks(log2([1,2,4])), xticklabels([1,2,4]), ylim([-.1 .4])
        title(locNames{iLoc})
    end
    
    set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)
    sgtitle(sprintf('linear slopes: %s - %d trials', subjName, nblocks/5*100), 'fontsize', 20)
end

