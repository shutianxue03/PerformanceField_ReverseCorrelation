
function error = fxn_getError_PR(iModelB, flagIncludePA, params, data)


% make prediction
pred = PR_pred(iModelB, params, data);
% fit trial-wise resonses
% only pC: from Fernandez2022
% plus pA: Dosher&Lu2008, Eq F2
pC = pred.metrics(3);
pA = pred.metrics(6);
m_pC = length(data.resp{1}) + length(data.resp{2}); % number of total trials
n_pC = sum(data.resp{1}) + sum(1-data.resp{2}); % number of total correct trials
m_pA = data.m_pA;% number of total pairs
n_pA = data.n_pA;% number of consistent pairs
fxn_nLL = @(p, m, n) -(n*log(p) + (m-n) * log(1-p));

if flagIncludePA, nLL = fxn_nLL(pC, m_pC, n_pC) + fxn_nLL(pA, m_pA, n_pA); 
else, nLL = fxn_nLL(pC, m_pC, n_pC); 
end

if nLL == Inf, nLL = 1/eps; end
error = nLL;
