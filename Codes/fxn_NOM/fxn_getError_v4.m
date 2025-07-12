
function nLL = fxn_getError_v4(iModelB, params, data)

% make prediction
% pred = PR_pred(iModelB, params, data); 
pred = PR_pred_v4(iModelB, params, data);  % predict pC and pA based on noisy IV with internal noises added on each trial

%% way 2: fit trial-wise responses
% only pC: from Fernandez2022
% plus pA: Dosher&Lu2008, Eq F2 in Appendix F (p39)
fxn_nLL = @(p, m, n) -(n*log(p) + (m-n) * log(1-p));% + sum(log(1:m)) - (sum(log(1:n)) - sum(log(1:(m - n))));
fxn_factorialLog =@(m,n) sum(log(1:m)) - sum(log(1:n)) - sum(log(1:(m-n)));

pC = pred.metrics(3);
m_pC = length(data.resp{1}) + length(data.resp{2}); % number of total trials
n_pC = sum(data.resp{1}) + sum(1-data.resp{2}); % number of total correct trials

pA = pred.metrics(6);
m_pA = length(data.respC(:, 1)); assert(m_pA == m_pC/2)% number of total pairs
n_pA = nansum(data.respC(:, 1));% number of consistent pairs

%%%%%%%%%%%%
% only fitting the absent trials
% pC = 1-pred.metrics(5); % pCR
% m_pC = length(data.resp{2}); % number of total trials
% n_pC = sum(1-data.resp{2}); % number of total correct trials
% 
% pA = pred.metrics(8);
% m_pA = length(data.respC(~isnan(data.respC(:, 3)), 3)); %assert(m_pA == m_pC/2)% number of total pairs
% n_pA = nansum(data.respC(:, 3));% number of consistent pairs
%%%%%%%

nLL = fxn_nLL(pC, m_pC, n_pC) + fxn_factorialLog(m_pC, n_pC) + fxn_nLL(pA, m_pA, n_pA) + fxn_factorialLog(m_pA, n_pA);

if nLL == Inf, nLL = 1/eps; end
