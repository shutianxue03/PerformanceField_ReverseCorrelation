

nbins = 15; % number of bars in the histogram
alpha = .3;

for iLoc = 1:nLoc
    figure('position', [0 200 800 nfilter_s*200])
    for itype = 2%1:ntypes
        switch itype
            case 1, istart = 1; iend = ntrials; % prs
            case 2, istart = ntrials+1; iend = ntrials*2; % abs
            case 3, istart = 1; iend = ntrials*2; % both
        end
        ntrials_ = length(istart: iend);
        x = energy_norm_allLoc{iLoc}(istart: iend, :);
        y = data_both(istart: iend, iLoc, 2);
        
        for ifilter_s = 1:nfilter_s
            ifilter_true = filter_s(ifilter_s);
            x_s = x(:, ifilter_true);
            %             fprintf('mean = %.4f, var = %.4f\n', mean(x_s), var(x_s))
            subplot(nfilter_s, ntypes, ntypes*(ifilter_s-1)+itype), hold on
            
            % plot the histogram of saying yes and saying no
            yyaxis left
            histogram(x_s(y==1), nbins, 'FaceColor', 'r', 'FaceAlpha', alpha, 'Normalization', 'probability', 'EdgeColor', 'none')
            histogram(x_s(y==0), nbins, 'FaceColor', 'b', 'FaceAlpha', alpha, 'Normalization', 'probability', 'EdgeColor', 'none')
            xline(mean(x_s(y==1)), 'r');
            xline(mean(x_s(y==0)), 'b');
            xline(mean(x_s), 'color', [1,0,1]);
            xline(0, 'k--');
            if ifilter_s+itype==2, ylabel('probability'), yticks([0,.25]), else, yticklabels([]), end
            set(gca,'ycolor',[.5,0, .5], 'ylim', [0,.25]);
            
            % plot the fitted logstist regression line
            yyaxis right
            plot(x_s, squeeze(yfit_all(itype, iLoc, ifilter_true, 1:ntrials_)), 'k.')
            if ifilter_s + itype == 2, ylabel('p(response = yes)'), yticks([0,.5, 1]), else, yticklabels([]),  end
            set(gca,'ycolor',[0,0,0], 'ylim', [0,1]);
            
            if ifilter_s + itype == 2, legend({'resp yes', 'resp no'}), xlabel('normalized energy'), xticks(-2:2:4)
            else, xticklabels([])
            end
            
            if ifilter_s== 1, title_type = sprintf('%s\n', typeNames{itype});
            else, title_type='';
            end
            xlim([-2,4])
            title(sprintf('%sfilter = %.2f, kernel = %.2f', title_type, filterSF_all(ifilter_true), kernel_allLoc(itype, iLoc, ifilter_true)))
        end
        
    end
    set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)
    sgtitle(locNames{iLoc}, 'fontsize', 25)
end

%%
figure('Position', [0 200 1200 300])
lineCol = {'r', 'b', 'k'};
for iLoc = 1:nLoc
    subplot(1, nLoc, iLoc), hold on
    for itype = 1:ntypes % prs, abs, both
        switch itype
            case 1, istart = 1; iend = ntrials;
            case 2, istart = ntrials+1; iend = ntrials*2;
            case 3, istart = 1; iend = ntrials*2;
        end
        x = energy_norm_allLoc{iLoc}(istart: iend, :);
        y = data_both(istart: iend, iLoc, 2);
        
        plot(filterSF_all_log, mean(x(y==1, :))-mean(x(y==0, :)), lineCol{itype})
        
    end
    yline(0); xline(log2(2));
    xticks(log2([1,2,4])), xticklabels([1,2,4])
    ylim([-.5,1])
    if iLoc==1,xlabel('SF channel (cpd)'), ylabel('mean energy (yes-no)'), end
    title(locNames{iLoc})
    
    set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)
end
