
% plot the psychometric function: p(YES) as a fxn of internal variable (IV)

figure('Position', [0 0 1000 500])
for itype = 1:2
    if itype == 1
        x_IV = IV_PRS_raw;
        y_resp = resp_test(boolean(iPRS_test));
    else
        x_IV = IV_ABS_raw;
        y_resp = resp_test(boolean(1-iPRS_test));
    end
    
    subplot(1,2,itype), hold on
    for step = [.01, .05, .1]
        quantiles = 0:step:1;
        x_bin = quantile(x_IV, quantiles);
        
        x=[]; y=[]; n=[];
        for ib = 1:length(quantiles)-1
            ind = x_IV>=x_bin(ib) & x_IV<=x_bin(ib+1);
            n(ib) = sum(ind);
            y_resp_ = y_resp(ind);
            y(ib) = mean(y_resp_);
            x_IV_ = x_IV(ind);
            x(ib) = mean(x_IV_);
        end
        
        % figure, plot(x_IV, y_resp, 'o')
        if step == .01
            plot(x, y, '+', 'markersize', 5)
        else
            plot(x, y, 'o', 'markersize', step*150)
        end
    end
    xlabel('IV')
    ylabel('P(YES)')
    ylim([0, 1])
    yticks(0:.1:1)
    %     xlim([min(x_IV), max(x_IV)])
    title(sprintf('Gabor-%s (n=%d)', namesType{itype}, length(x_IV)))
    grid on
end
sgtitle(sprintf('%s-%s\n[A%d] %s', subjName, namesLocComb{iLocComb}, iModelA, namesModelA{iModelA}))
set(findall(gcf, '-property', 'fontsize'), 'fontsize',20)
set(findall(gcf, '-property', 'linewidth'), 'linewidth',1.5)

