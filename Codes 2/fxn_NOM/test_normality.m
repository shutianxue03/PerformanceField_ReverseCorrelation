% function test_normality

% conclusion: most linear/log IV distributions are skewed...
clc
close all
clear all

isubj = 1;
subjList = {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'AS', 'HL', 'FH', 'HA', 'CS', 'DT'};
nblocks_allSubj = [200, 240, 210, 240, 220, 220, 215, 220, 205, 210, 205, 205];

nB = 1;
nfiltersOri = 29;
nfiltersSF = 29;
itype = 3;
iLoc = 1;
ORI_bound = [7, 23];

ind_all = combvec(2:3, 1:2, 1:2, 0:1);
nind = size(ind_all, 2);

titles{1} = {'rawT', 'posT', 'reconT'};
titles{2} = {'dotP', 'Conv'};
titles{3} = {'sum', 'max'};
titles{4} = {'Linear', 'Log'};

subjName_ = subjList{isubj};
nblocks = nblocks_allSubj(isubj);

fprintf([subjName_, ' Loading...'])
load(sprintf('Data_OOD/%s%d_B%d_%d_%d.mat', subjName_, nblocks, nB, nfiltersOri, nfiltersSF))
load(sprintf('Data_OOD/%s%d_energy_%d_%d.mat', subjName_, nblocks, nfiltersOri, nfiltersSF))
load(sprintf('Data_OOD/%s%d_behavMeas.mat', subjName_, nblocks))
fprintf('DONE\n')

%%
for iind = 1:nind
    templateType.i = ind_all(1, iind); title1 = titles{1}{templateType.i};
    convolveType.i = ind_all(2, iind);title2 = titles{2}{convolveType.i};
    IVType.i = ind_all(3, iind);title3 = titles{3}{IVType.i};
    logFlag = ind_all(4, iind);title4 = titles{4}{logFlag+1};
    
    %% get template and energy of that location
    [template, e2D] = fxn_getTemplate(kernels2D_perComb, itype, iLoc, templateType, energy2D_allT_perComb);
    
    %% get IV
    [IV_PRS_, IV_ABS_, imax_allT] = fxn_getIV(ntrials, e2D, convolveType, IVType, template, ORI_bound);
    if sum(IV_PRS_>0)<length(IV_PRS_), fprintf('[PRS] There are negative IV\n'), end
    if sum(IV_ABS_>0)<length(IV_ABS_), fprintf('[ABS] There are negative IV\n'), end

    %% test if IV is normally distributed (Anderson-Darling test)
    IV_PRS_ln = IV_PRS_;
    IV_ABS_ln = IV_ABS_;
    IV_PRS_log = log2(IV_PRS_); 
    IV_ABS_log = log2(IV_ABS_); 
    
    if logFlag
        IV_PRS = IV_PRS_log;
        IV_ABS = IV_ABS_log;
    else
        IV_PRS = IV_PRS_ln;
        IV_ABS = IV_ABS_ln;
    end
    [hPRS, pPRS] = swtest(IV_PRS, .2); if hPRS == 0, sPRS = 'succeeded'; else, sPRS = 'FAILED'; end
    [hABS, pABS] = swtest(IV_ABS); if hABS == 0, sABS = 'succeeded'; else, sABS = 'FAILED'; end
    
    %% plot IV distribution
    figure('Position', [2000 0 500 500])
    subplot(2,2,[1,2]), hold on
    histogram(IV_PRS, 50, 'DisplayStyle', 'stairs', 'EdgeColor', 'r')
    xline(median(IV_PRS), 'r-', 'linewidth', 2);
    xline(mean(IV_PRS), 'r--', 'linewidth', 2);
    
    histogram(IV_ABS, 50, 'DisplayStyle', 'stairs', 'EdgeColor', 'b')
    xline(median(IV_ABS), 'b-', 'linewidth', 2);
    xline(mean(IV_ABS), 'b--', 'linewidth', 2);
    legend({'PRS', 'med', 'mean', 'ABS', 'med', 'mean'})
    
    subplot(2,2,3), qqplot(IV_PRS), title(sprintf('[PRS] Normality test %s (p = %.3f)', sPRS, pPRS))
    subplot(2,2,4), qqplot(IV_ABS), title(sprintf('[ABS] Normality test %s (p = %.3f)', sABS, pABS))
    sgtitle(sprintf('%s\n%s-%s-%s-%s\n[PRS] med = %.2f, SD = %.3f // [ABS] med = %.2f, SD = %.3f', subjName, title1, title2, title3, title4, median(IV_PRS), std(IV_PRS), median(IV_ABS), std(IV_ABS)))
end % end of iind

