
%% instruction
% make sure someone's data is hangingin there
% set a breakpoint at 'for ifilterOri = 1:nfiltersOri' in SX_sim07_RC.m
% Then run SX_RC6_RC
% Then copy and paste this script in the command window 
% (make sure the Link function is logit, not probit, since JASP runs logistic regression)
% Open JASP and open file 'test.csv'

ifilterOri = 10;   
ifilterSF= 8;  
x = energy_both(:, ifilterOri, ifilterSF); 
csvwrite('Data/test.csv', [[1;x],[1;y]]);
[beta,~,stats] = glmfit(x,y,'binomial','Link','logit');
R2 = 1-var(stats.resid)./var(y); 
pValues = stats.p;
yfit = glmval(beta,x,'logit');

% get R2_Tjur
SS_res = sumsqr(stats.resid);
SS_mod = sumsqr(mean(y) - yfit);
SS_total = sumsqr(mean(y) - y); % unquivalent to SS_res+SS_mod
R2_mod = SS_mod/SS_total;
R2_res = 1-SS_res/SS_total;
R2_Tjur = (R2_mod+R2_res)/2;

fprintf('R2 = %.3f, R2 Tjur = %.3f, intercept = %.3f, slope = %.3f\n', R2, R2_Tjur, beta(1), beta(2))