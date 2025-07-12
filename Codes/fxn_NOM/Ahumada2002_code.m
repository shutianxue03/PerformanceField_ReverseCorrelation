

pCS02 = [911  417  149    180  352   941] ;  
% first 3: given a signal-abs stim; last 3: given a signal-prs stim
% for each group: both respond abs, one abs and one prs, both respond prs

pCS16 = [1290  552  115    370  657  1016] ; 

pCS = pCS16;
pFA = (pCS(3) + pCS(2)/2)/sum(pCS(1:3)); % 0.2420
pHit = (pCS(6) + pCS(5)/2)/sum(pCS(4:6)); % 0.7583

dH = norminv(pHit) - norminv(pFA); % 1.4006
d0 = 3.2006; % stim SNR % how to determine this?
alpha2_dprime = (dH/d0)^2; % 0.1915, alpha^2, eq 3.2.1 Ahumada2002 and eq8 in Burgess&Colborne1987
% criterion = -norminv(pFA); 
criterion = -(norminv(pHit)+norminv(pFA));
nPairs =  [sum(pCS(1:3)), sum(pCS(4:6))];
pA_both = [(pCS(1)+pCS(3))/nPairs(1), (pCS(4)+pCS(6))/nPairs(2)]; % response consistency?

pA_CI = nan(2,3); % confidence interval
for s = 1:2 % tgt-prs and abs
    pA = pA_both(s);
    n = nPairs(s);
    CI_step = 1.96*sqrt(((1-pA)*pA)/n); % 1.96: 95%
    pA_CI(s,:) = [pA, pA - CI_step, pA + CI_step]; % mid, lb, ub of the CI
end

%% setimate noise SD ratio, i.e., alpha^2
n = 25; % integrtion approx. parameter
range = 4; % ditto

lb = .001; % a minimal number to avoid using 0
ub = 1;
mu_both = [0, dH]; % signal-abs, signal-prs

sigmaI_est = nan(2,3);
for s = 1:2
    mu = mu_both(s);
    for i = 1:3
        pA = pA_CI(s,i); % the measured P(A)
        % to estimate internal/external ratio that predicts psame
        f = @(sigma_I) (pA-probSame(sigma_I, criterion, mu, n, range))^2; 
        sigmaI_est(s,i) = fminbnd(f, lb, ub);
    end
end

alpha2_est = 1-sigmaI_est.*sigmaI_est;

%%
figure, hold on
errorbar(1,alpha2_est(1,1),alpha2_est(1,1)-alpha2_est(1,2), 'ko')
errorbar(3,alpha2_est(2,1),alpha2_est(2,1)-alpha2_est(2,2), 'ko')
plot(2, alpha2_dprime,'kx')

xlim([.5, 3.5])
xticks(1:3), xticklabels({'stim-abs', 'd', 'stim-prs'})
grid on
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(findall(gcf,'-property','linewidth'),'linewidth',1.5)


%%
function prob = probSame(sigma_I, criterion, mu, n, range)
% seems to be eq5 in Burgess & Colborne (1988)
% seems to be eq xx in Richards & Zhu (1994)
sigma_E = sqrt(1-sigma_I^2); % assuming sigma_E^2+sigma_I^2 = 1 (?)
a = (mu-criterion)/sigma_I;
k = sigma_E/sigma_I; % Burgess & Colborne, 1988, eq5, p4, top
t = -range : 1/n : range;
x = a + k * t; % input of fxn Q(), Burgess & Colborne, 1988
p = normcdf(x);
p = p .* (1-p) .* exp(-t.^2/2);
prob = 1-2*sum(p)/(n*sqrt(2*pi));

end

    

