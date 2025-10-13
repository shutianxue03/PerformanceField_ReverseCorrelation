
%%
% compare metrics (dprime, crterion, pC3, pA3) calculated by 
% 1. math 
% 2. simulation

% the simulated pA is very different from the mathematically predicted pA
% thus I decided to simulate pA when fitting

%%
clc
close all
% predict pA from simulation and math
thresh = .25;
mu_IV_PRS = .31;
mu_IV_ABS = .23;
std_IV_PRS = .05;
std_IV_ABS = .04;

ntrials_sim = 1e6;
ni = 100;

%% plot the theoretical IV distribution
xmin = 0;
xmax = .5;
x = linspace(xmin, xmax, 1e5);
y_PRS = normpdf(x, mu_IV_PRS, std_IV_PRS);
y_ABS = normpdf(x, mu_IV_ABS, std_IV_ABS);

% plot the IV distribution
figure('Position', [200 200 800 300])
subplot(1,2,1), hold on
plot(x, y_PRS, 'r-', 'linewidth', 1.5)
plot(x, y_ABS, 'b-', 'linewidth', 1.5)
xline(mu_IV_PRS, 'r', 'linewidth', 1.5);
xline(mu_IV_ABS, 'b', 'linewidth', 1.5);
xline(thresh, 'k', 'linewidth', 2);
errorbar(mu_IV_PRS, max(y_PRS)/2, std_IV_PRS, 'r-', 'horizontal', 'linewidth', 1.5)
errorbar(mu_IV_ABS, max(y_ABS)/2, std_IV_ABS, 'b-', 'horizontal', 'linewidth', 1.5)

%% predict dprime/criterion/pHit/pFA/pC by math
pHit_math = 1-normcdf(thresh, mu_IV_PRS, std_IV_PRS);
pFA_math = 1-normcdf(thresh, mu_IV_ABS, std_IV_ABS);
pC_math = (pHit_math+1-pFA_math)/2;
[dprime_math, criterion_math] = SX_sim06_SDT(pHit_math, pFA_math, 1e5); % 1e5 is just an arb number to ensure d'/c not to be Inf

%% predict pA by math
% SX_normPDF = 
pdf_PRS_ = normpdf(x-thresh, mu_IV_PRS, std_IV_PRS); pdf_PRS = pdf_PRS_/sum(pdf_PRS_);
pdf_ABS_ = normpdf(x-thresh, mu_IV_ABS, std_IV_ABS);pdf_ABS = pdf_ABS_/sum(pdf_ABS_);
cdf_PRS = normcdf(x-thresh, mu_IV_PRS, std_IV_PRS);
cdf_ABS = normcdf(x-thresh, mu_IV_PRS, std_IV_PRS);


pA_PRS_math = sum((cdf_PRS.^2 + (1-cdf_PRS).^2) .* pdf_PRS);
pA_ABS_math =  sum(((1-cdf_ABS).^2 + cdf_ABS.^2) .* pdf_ABS);
pA_math = mean([pA_PRS_math, pA_ABS_math]);

pHit_math = sum((1-cdf_PRS)  .* pdf_PRS);
pFA_math = sum((1-cdf_ABS) .* pdf_ABS);
pC_math = (pHit_math+1-pFA_math)/2;

%% predict all metrics by simulating data
% simulate the internal responses (IR) for both passes
dprime_sim_allB = nan(ni, 1);
criterion_sim_allB= dprime_sim_allB;
pC_sim_allB = dprime_sim_allB;
pHit_sim_allB = dprime_sim_allB;
pFA_sim_allB = dprime_sim_allB;
pA_PRS_sim_allB = dprime_sim_allB;
pA_ABS_sim_allB= dprime_sim_allB;
pA_sim_allB =  dprime_sim_allB;

parfor ii = 1:ni
    % generate IV of multiple trials
    IV_PRS_ = randn(2, ntrials_sim) * std_IV_PRS + mu_IV_PRS;
    IV_ABS_ = randn(2, ntrials_sim) * std_IV_ABS + mu_IV_ABS;
   
    % get response (1=YES , 0=NO)
    resp_PRS = IV_PRS_ >= thresh;
    resp_ABS = IV_ABS_ >= thresh;
    
    % calculate dprime/criterion/pC3
    pHit_sim = mean(resp_PRS(:));
    pFA_sim = mean(resp_ABS(:)); 
    pC_sim = mean([pHit_sim, 1-pFA_sim]);
    [dprime_sim, criterion_sim] =  SX_sim06_SDT(pHit_sim, pFA_sim);
    
    % calculate pA3
    pA_PRS_sim = mean(resp_PRS(1,:) == resp_PRS(2,:));
    pA_ABS_sim = mean(resp_ABS(1,:) == resp_ABS(2,:));
    pA_sim = mean([resp_PRS(1,:) == resp_PRS(2,:), resp_ABS(1,:) == resp_ABS(2,:)]);
    
    dprime_sim_allB(ii) = dprime_sim;
    criterion_sim_allB(ii) = criterion_sim;
        
    pC_sim_allB(ii) = pC_sim;
    pHit_sim_allB(ii) = pHit_sim;
    pFA_sim_allB(ii) = pFA_sim;
    
    pA_sim_allB(ii) = pA_sim;
    pA_PRS_sim_allB(ii) = pA_PRS_sim;
    pA_ABS_sim_allB(ii) = pA_ABS_sim;

end

%% plot one simulated IV sample
IV_PRS = randn(2, ntrials_sim) * std_IV_PRS + mu_IV_PRS; % IR_PRS is created inside the parfor thus is not available outside the loop
IV_ABS = randn(2, ntrials_sim) * std_IV_ABS + mu_IV_ABS;

subplot(1,2,2), hold on,
histogram(IV_PRS, 'FaceColor', 'r', 'FaceAlpha', .3, 'EdgeColor', 'none', 'Normalization', 'probability');
histogram(IV_ABS, 'FaceColor', 'b', 'FaceAlpha', .3, 'EdgeColor', 'none', 'Normalization', 'probability');
xline(thresh, 'linewidth', 2);
xlim([xmin, xmax])
title(sprintf('sim dprime = %.4f (%.4f)\nsim criterion = %.4f (%.4f)', mean(dprime_sim_allB), std(dprime_sim_allB), mean(criterion_sim_allB), std(criterion_sim_allB)))

%%
figure('Position', [0 200 1200 200])
plot_compare(dprime_math, dprime_sim_allB, 1, 'dprime')
plot_compare(criterion_math, criterion_sim_allB, 2, 'criterion')
plot_compare(pC_math, pC_sim_allB, 3, 'pC')
plot_compare(pHit_math, pHit_sim_allB, 4, 'pHit')
plot_compare(pFA_math, pFA_sim_allB, 5, 'pFA')
plot_compare(pA_math, pA_sim_allB, 6, 'pA')
plot_compare(pA_PRS_math, pA_PRS_sim_allB, 7, 'pA PRS')
plot_compare(pA_ABS_math, pA_ABS_sim_allB, 8, 'pA ABS')


%%
function plot_compare(data_math, data_sim_allB, iplot, title_)
data_sim_mean = mean(data_sim_allB);
data_sim_med = median(data_sim_allB);
% data_sim_std = std(data_sim_allB);
% data_sim_CI_neg = data_sim_med - quantile(data_sim_allB, .16);
% data_sim_CI_pos = quantile(data_sim_allB, .84) - data_sim_med;

subplot(1,8,iplot), hold on
xline(data_math, 'k', 'LineWidth', 2);
histogram(data_sim_allB, 'FaceColor', 'k', 'FaceAlpha', .3, 'EdgeColor', 'none', 'Normalization', 'probability')
% bar([1,2,3], [data_math, data_sim_mean, data_sim_med])
xline(data_sim_mean, 'r-');
% errorbar(data_sim_mean, .2, data_sim_std, 'k', 'horizontal', 'linewidth', 1, 'CapSize', 0)
xline(data_sim_med, 'b-');
% errorbar(data_sim_med, .1, data_sim_CI_neg, data_sim_CI_pos, 'k', 'horizontal','linewidth', 1, 'CapSize', 0)
% xticks(1:3)
% xticklabels({'MATH', 'SIM (mean)', 'SIM (median)'})
title(title_)
% ylim([0, 1])
% yline(.5);
if iplot == 1, legend({'Math', 'Simulation', 'Mean', 'Median'}, 'Location', 'south'),end
% set(findall(gcf, '-property', 'FontSize'), 'FontSize',18)
% sgtitle(sprintf('PRS ~ N(%.1f, %.1f)\nABS ~ N(%.1f, %.1f)\nthresh = %.1f\n #sim trials = 1e%d', mu_IV_PRS, std_IV_PRS, mu_IV_ABS, std_IV_ABS, thresh, log10(ntrials_sim)))

end