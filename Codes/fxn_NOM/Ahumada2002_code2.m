

f = @(gamma, p, betaH, dprime_e, elim) (p - probR00(gamma, betaH, dprime_e, elim))^2;

dprime_e = 0.04; % why 0.04?
elim = 4.0; 
gamma_lb = 0.001; 
gamma_ub = 100; 
nP = [646, 266, 89; 137, 218, 644]; 
nP0 = sum(nP(1,:)); % number of signal-abs trials; 
nP1 = sum(nP(2,:)); % number of signal-prs trials; 
P000 = nP(1,1)/nP0; % prop of both resp = abs given signal abs; 
P100 = nP(2,1)/nP1;  % prop of both resp = abs given signal prs; 
P00 = (nP(1,1) + nP(1,2)/2)/nP0; % prop of resp = abs given signal abs; 
P10 = (nP(2,1) + nP(2,2)/2)/nP1;  % prop of resp = abs given signal prs; 

% independence predictions (alpha = 0) 
p000_pred=P00*P00;
p100_pred=P10*P10; 
pCI(1,:) = [P000, confIntP(P000,nP0,1.96)]; 
pCI(2,:) = [P100, confIntP(P100,nP1,1.96)]; 

betaH(1) = norminv(P00); 
betaH(2) = norminv(P10);

gamma_est = nan(2,3);
for i=1:2 % abs and prs
    for j=1:3 % CI
        gamma_est(i,j) = fminbnd(f, gamma_lb, gamma_ub, [], pCI(i,j), betaH(i), dprime_e, elim); 
    end
end

alpha2 = 1 ./ (1 + gamma_est .* gamma_est);

%% subfunctions
function prob = probR00(gamma,betaH,dprime_e,elim) 
% Prob{both responses are 0} from gamma, betaH 
e =-elim:dprime_e:elim; 
p = normcdf((betaH * sqrt(1+gamma^2) - e)/gamma); % Ahumada (2002) eq3.1.4
p = p.^2.*exp(e^2/2); 
prob = sum(p)*(dprime_e/sqrt(2*pi));  % why this??
end

function prob = probR11(gamma,betaH,dprime_e,elim) 
%Prob{both responses are 1} from gamma, betaH 
%Ahumada (2002) Equation (3.1.5) 
e =-elim:dprime_e:elim; 
p = 1-normcdf((betaH * sqrt(1+gamma^2) - e)/gamma); % Ahumada (2002) eq3.1.5
p = p.^2.*exp(e^2/2); 
prob = sum(p)*(dprime_e/sqrt(2*pi));  % why this??
end

function interval = confIntP(p, n, z) 
% confidence interval for a proportion
% p1+z*sqrt(p1*(1-p1)/n) = p 
%(z^2/n)*(p1-p1^2) = (p-p1)^2 
%(z^2/n+1)*p1^2 +(-2*p-z^2/n)*p1 +(p^2) = 0 
% x = (b +- sqrt(b^2-4*a*c))/(2*a) 

b = 2*p+z^2/n; a = z^2/n+1 ; 
c = p*p; 
c = sqrt(b^2-4*a*c)/(2*a); 
b = b/(2*a); 
interval = [b - c , b+c]; 
%p = .1; n = 100 ;z = 1.96;
end


