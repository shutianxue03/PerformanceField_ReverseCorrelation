
function error = fxn_getError_v3(iModelB, params, data)

% make prediction
pred = PR_pred_v3(iModelB, params, data);

% calculate nLL
pC = pred.metrics(3);
m = length(data.resp{1}) + length(data.resp{2}); % number of total trials
n_PRS = sum(data.resp{1}); % number of correct PRS trials
n_ABS = length((data.resp{2})) - sum(data.resp{2}); % number of correct ABS trials
n = n_PRS + n_ABS; % number of total correct trials
fxn_nLL = @(pC, m, n) -(n*log(pC) + (m-n) * log(1-pC));
nLL = fxn_nLL(pC, m, n); if nLL == Inf, nLL = 1/eps; end
error = nLL;
