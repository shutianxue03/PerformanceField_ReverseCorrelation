

PF_titles = {'LR diff', 'HVA', 'VMA'};
nPFtitles = length(PF_titles);
xlabels = {'Left', 'VM', 'Upper VM'};
ylabels = {'Right', 'HM', 'Lower VM'};
ncomp = length(xlabels);
SDT_PFticks = {[.5,1,1.5]; [-1,0,1]; [0, .2, .4]};
ori_PFticks = {5:5:25, 0:1:6, [-.025, .025]};
SF_PFticks = {[1.5, 2.5], [1,4], [0,4], [-4,0]};

figure('Position', [0 0 800 900])
% just combine 3 SDT data together
dataSDT = cat(3, dprime_allSubj, criterion_allSubj, exp(RT_log_allSubj));

itype = 3;
dataTypes = {'SDT'};
%     dataTypes = {'SDT', 'Ori params', 'SF params'};
for dd = 1:length(dataTypes)
    dataType = dataTypes{dd};
    switch dd
        case 1, measNames = namesSDT; axisticks = SDT_PFticks;
        case 2, measNames = modelParamsNamesOri; axisticks = ori_PFticks;
        case 3, measNames = modelParamsNamesSF; axisticks = SF_PFticks;
    end
    
    nMeas = size(dataSDT, 3);
    iplotInd = reshape(1:nPFtitles * nMeas, nMeas, nPFtitles).';
    

    % data: nsubj x nLoc x nparams
    switch dd
        case 1, data = dataSDT; axisticks = SDT_PFticks;
        case 2, data = squeeze(paramsOri_allSubj(:, itype, :, :)); axisticks = ori_PFticks;
        case 3, data = squeeze(paramsSF_allSubj(:, itype, :, :)); axisticks = SF_PFticks;
    end
    
    for n = 1:nMeas % nMeas: measurement of each data set ()
        
        Left = squeeze(data(:,2,n));
        Right = squeeze(data(:,4,n));
        UVM = squeeze(data(:,3,n));
        LVM = squeeze(data(:,5,n));
        HM = nanmean([Left, Right],2);
        VM = nanmean([UVM, LVM],2);
        
        for ip = 1:ncomp
            
            switch ip, case 1, x = Left; y = Right; case 2, x = VM; y = HM; case 3, x = LVM; y = UVM;end
            x_mean = mean(x);
            y_mean = mean(y);
            
            subplot(ncomp , nMeas, iplotInd(ip, n)), hold on
            
            % plot idvd data
            for isubj_p = 1:nsubj, plot(x(isubj_p),y(isubj_p), [marks_allSubj{isubj_p}, 'k']), end
            
            axisticks_ = axisticks{n};
            axismin = axisticks_(1);
            axismax = axisticks_(end);
            
            % plot mean and error
            plot(x_mean, y_mean, 'k.') % mean
            errorbar(x_mean, y_mean, std(x)/sqrt(nsubj), 'horizontal', 'k')
            errorbar(x_mean, y_mean, std(y)/sqrt(nsubj), 'vertical', 'k')
            
            % plot the diagonal line
            plot([axismin, axismax], [axismin, axismax], '-', 'color', [.5,.5,.5]) % the diagonal line
            xlim([axismin, axismax]), ylim([axismin, axismax])
            xticks(axisticks_)
            yticks(axisticks_)
            axis square
            box on
            
            if ip == 1 , title(measNames{n}, 'FontSize',15), end
            if n == 1
                xlabel(xlabels{ip}, 'FontSize',15), ylabel(ylabels{ip}, 'FontSize',15)
                xline(dprime_theo, 'color', [.5,.5,.5]); yline(dprime_theo, 'color', [.5,.5,.5]);
            elseif n ==2, xline(0, 'color', [.5,.5,.5]); yline(0, 'color', [.5,.5,.5]);
            end
            if ip + n == 2, legend(subjList, 'Location', 'best'), end
        end
    end
    
    set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)
    sgtitle(sprintf('n=%d, %s', nsubj, dataTypes{dd}), 'FontSize',25)
    
end