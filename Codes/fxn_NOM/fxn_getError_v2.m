
function nLL = fxn_getError_v2(iModelB, flagIncludePA, params, data)


% make prediction
pred = PR_pred(iModelB, flagIncludePA, params, data);

%% way 1: get nLL for PRS and ABS separately
% nLL = nan(1,2);
% for i = 1:2
%     if i==1 % PRS
%         pC = pred.metrics(4);
%         pA = pred.metrics(7);
%     else
%         pC = 1-pred.metrics(5);
%         pA = pred.metrics(8);
%     end
%     m_pC = length(data.resp{i}) ;
%     n_pC = sum(data.resp{i});
%     m_pA = length(data.respC(:, i+1));% number of total pairs
%     n_pA = nansum(data.respC(:, i+1));% number of consistent pairs
%     fxn_nLL = @(pC, m, n) -(n*log(pC) + (m-n) * log(1-pC));% + sum(log(1:m)) - (sum(log(1:n)) - sum(log(1:(m - n))));
%     nLL(i) = fxn_nLL(pC, m_pC, n_pC);% + fxn_nLL(pA, m_pA, n_pA); if nLL == Inf, nLL = 1/eps; end
% end
% nLL = sum(nLL);

%% way 2: fit trial-wise responses
% only pC: from Fernandez2022
% plus pA: Dosher&Lu2008, Eq F2
pC = pred.metrics(3);
pA = pred.metrics(6);
% [pC, pA] = getpC_pA(data.cst, .2, sigma_int);

m_pC = length(data.resp{1}) + length(data.resp{2}); % number of total trials
n_pC = sum(data.resp{1}) + sum(1-data.resp{2}); % number of total correct trials
m_pA = length(data.respC(:, 1)); assert(m_pA == m_pC/2)% number of total pairs
n_pA = nansum(data.respC(:, 1));% number of consistent pairs
fxn_nLL = @(p, m, n) -(n*log(p) + (m-n) * log(1-p));% + sum(log(1:m)) - (sum(log(1:n)) - sum(log(1:(m - n))));
if flagIncludePA, nLL = fxn_nLL(pC, m_pC, n_pC) + fxn_nLL(pA, m_pA, n_pA);
else, nLL = fxn_nLL(pC, m_pC, n_pC);
end
if nLL == Inf, nLL = 1/eps; end

%% helper fxn - getpC_pA
    function [pC, pA] = getpC_pA(cst, sigma_ext, sigma_int)
        % SX_normPDF = @(x,mu,sigma) exp(-(x-mu).^2/(2*sigma^2));
        
        % get x
        x = (-1e4:1e4)/100;
        sigma_int_PRS = sigma_int;
        sigma_int_ABS = sigma_int;
        
        % combine internal & external noise
        sigma_IE = sqrt(sigma_ext^2+(sigma_int_PRS^2+sigma_int_ABS^2)/2);
        t = cst + x * sigma_IE;
        
        % get pC
        pC = sum((sigma_IE)*normpdf(t,cst,sigma_int_PRS).*normcdf(t,0,sigma_int_ABS)/100);
        
        sigma_ext2 = sigma_ext * sqrt(2);
        sigma_int2 = sqrt(sigma_int_PRS^2+sigma_int_ABS^2);
        t = cst + x * sigma_ext2;
        
        % get pA
        pA = sum((sigma_ext2)*normpdf(t,cst,sigma_ext2).*(normcdf(t,0,sigma_int2).^2+(1-normcdf(t,0,sigma_int2)).^2)/100);
        
    end
end
