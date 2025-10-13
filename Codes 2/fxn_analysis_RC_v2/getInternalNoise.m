
function [alpha2_est, alpha2_dprime] = getInternalNoise(pCS, d0)

assert(length(pCS) == 6)
% pCS:
% first 3: given a signal-abs stim; last 3: given a signal-prs stim
% for each group: both respond abs, one abs and one prs, both respond prs

pFA = (pCS(3) + pCS(2)/2)/sum(pCS(1:3)); % 0.2420
pHit = (pCS(6) + pCS(5)/2)/sum(pCS(4:6)); % 0.7583

dprimeH = norminv(pHit) - norminv(pFA); % 1.4006
alpha2_dprime = (dprimeH/d0)^2; % 0.1915, alpha^2, eq 3.2.1 Ahumada2002 and eq8 in Burgess&Colborne1987
beta0 = -norminv(pFA); % bias, eq2.2.8, Ahumada 2002
nPairs =  [sum(pCS(1:3)), sum(pCS(4:6))];
pA_both = [(pCS(1)+pCS(3))/nPairs(1), (pCS(4)+pCS(6))/nPairs(2)]; % response consistency for signal-abs and prs

% pA_CI = nan(2,3); % confidence interval
% for s = 1:2 % tgt-prs and abs
%     pA = pA_both(s);
%     n = nPairs(s);
%     CI_step = 1.96*sqrt(((1-pA)*pA)/n); % 1.96: 95% confidence interval
%     pA_CI(s,:) = [pA, pA - CI_step, pA + CI_step]; % mid, lb, ub of the CI
% end

%% etimate noise SD ratio, i.e., alpha^2
lb = .001; % a minimal number to avoid using 0
ub = 1;
mu_both = [0, dprimeH]; % signal-abs, signal-prs

s_est = nan(1,length(mu_both0));
for s = 1:2
    mu = mu_both(s);    
    pA = pA_both(s);
    fxn = @(sigma_i) (pA - probSame(sigma_i, beta0, mu))^2;
    s_est(s) = fminbnd(fxn, lb, ub);
end

alpha2_est = 1-s_est.*s_est;

%%
function prob = probSame(sigma_I, beta0, mu)
    % seems to be eq5 in Burgess & Colborne (1988)
    % seems to be eq xx in Richards & Zhu (1994)
    n = 25; % integrtion approx. parameter
    range = 4; % ditto

    sigma_E = sqrt(1-sigma_I^2); % assuming sigma_E^2+sigma_I^2 = 1 (?)
    a = (mu-beta0)/sigma_I;
    k = sigma_E/sigma_I; % Burgess & Colborne, 1988, eq5, p4, top
    t = -range : 1/n : range;
    x = a + k * t; % input of fxn Q(), Burgess & Colborne, 1988
    p = normcdf(x);
    p = p .* (1-p) .* exp(-t.^2/2);
    prob = 1-2*sum(p)/(n*sqrt(2*pi));
end


end

