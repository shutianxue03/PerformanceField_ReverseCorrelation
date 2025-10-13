%% get internal noise and efficiency
iCuedLoc_last5 = iCuedLoc(end-4:end);
cst_allB_last5 = cst_perB(end-4:end);
InternalNoise_allLoc = nan(nLoc, 2); % signal-abs and -prs
InternalNoise_dprime_all = nan(nLoc, 1);
efficiency_allLoc = nan(nLoc, 2); % signal-abs and -prs

sigma_ext = 1;
sigma0 = .3; % initial value
sigma_lb = .1;
sigma_ub = 1;
sigmaINT_est_allLoc = nan(1,nLoc);
options = optimoptions('fmincon','MaxIterations',5000,'Display','off');

%%
figure('Position', [0 200 1500 250]) % plot data (averaged over signal cst) and the estimated pA and pC
for iLoc = 1:nLoc
    
    signalCST_all = cst_unik_perLoc{iLoc}.';
    nsignalCST_all = length(signalCST_all);
    pA_data = pA_unik_perLoc{iLoc}; % pA calculated based on the unique contrast at each loc
    pC_data = pC_unik_perLoc{iLoc};
    fxn = @(sigmaINT) getDiff(signalCST_all, sigma_ext, sigmaINT, pA_data, pC_data);
    sigmaINT_est = fmincon(fxn, sigma0, [], [], [], [], sigma_lb, sigma_ub, [], options);
    sigmaINT_est_allLoc(iLoc) = sigmaINT_est;
    
    % predict pC and pA
    pC_est_all = nan(1, nsignalCST_all);
    pA_est_all = pC_est_all;
    for iCST = 1:nsignalCST_all
        signalCST = signalCST_all(iCST);
        [pC_est, pA_est] = getpC_pA(signalCST, sigma_ext, sigmaINT_est);
        pC_est_all(iCST) = pC_est;
        pA_est_all(iCST) = pA_est;
    end
    
    % plot
    subplot(1,nLoc, iLoc), grid on, hold on, box on
    
    for iCST = 1:nsignalCST_all
        % averaged across cst
        markerSize = nblocks_perCST_perLoc{iLoc}(iCST)/sum(nblocks_perCST_perLoc{iLoc})*20+10;
        plot(pA_data(iCST), pC_data(iCST), 'o', 'color', colorUniks(iCST, :), 'markersize', markerSize)
    end
    plot(pA_est_all', pC_est_all, '.-k')
    xlim([.5, 1]), xticks(.5:.1:1), xlabel('pA')
    ylim([.5, 1]), yticks(.5:.1:1), ylabel('pC')
    title(sprintf('%s, \\sigma_{int}=%.4f', locNames{iLoc}, sigmaINT_est))
end
set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)
set(findall(gcf, '-property', 'linewidth'), 'linewidth',1.5)

%% plot the estimated internal noise SD across loc
polarAngFull = (0:90:360)/180*pi;
polarIndFull = [4,3,2,5,4];
figure('Position', [0 200 800 250])
subplot(1,2,1), hold on
bar(sigmaINT_est_allLoc, 'handlevisibility', 'off')
errorbar(1,mean(sigmaINT_est_allLoc(2:end)), std(sigmaINT_est_allLoc(2:end))/sqrt(4), 'o')
xticks(1:5), xticklabels(locNames), xtickangle(45)
ylabel('SD of internal noise')
legend('periphery ave')

subplot(1,2,2), hold on
polarAxesHandle = subplot(1,2,2);
ax = gca;
ax.XTick = [];
ax.YTick = [];
polaraxes('Units',polarAxesHandle.Units,'Position',polarAxesHandle.Position)
hold on
polarplot((0:90:360)/180*pi, sigmaINT_est_allLoc(polarIndFull)/sigmaINT_est_allLoc(1), 'k')
ax = gca;
ax.ThetaTick = [];
% rlim([.8, 1])
title('relative to the center')

set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)
set(findall(gcf, '-property', 'linewidth'), 'linewidth',1.5)

%% helper functions
function diff = getDiff(signalCST_all, sigma_ext, sigmaINT, pA_data, pC_data)
% calculate the squared error between data and estimated pA and pC
nsignalCST_all = length(signalCST_all);
diff = 0;
for iCST = 1:nsignalCST_all
    [pC_est, pA_est] = getpC_pA(signalCST_all(iCST), sigma_ext, sigmaINT);
    diff = diff + sum((pA_data - pA_est).^2) + sum((pC_data - pC_est).^2);
end
end

function [pC, pA] = getpC_pA(signalCST, sigmal_ext, sigmaINT)
% calculate pC and pA given signal contrast, external noise SD and internal
% noise SD
sigma_abs = sigmaINT; % sigmaINT is sigma_prs
x_range = (-1e4:1e4)/100;

N = sqrt(sigmal_ext^2+(sigmaINT^2 + sigma_abs^2)/2);
t = signalCST + x_range * N;
pC = sum(N/100*normpdf(t,signalCST,sigmaINT).*normcdf(t,0,sigma_abs));

sX = sigmal_ext * sqrt(2);
b = sqrt(sigmaINT ^ 2 + sigma_abs ^ 2);
t = signalCST + x_range * sX;
pA = sum(sX/100*normpdf(t,signalCST,sX).*(normcdf(t,0,b).^2+(1-normcdf(t,0,b)).^2));

end