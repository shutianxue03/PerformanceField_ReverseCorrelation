
function [ub, lb, params0] = MC_getLimits(paramInd, ub_full, lb_full, nLocTrue)
ub = [];
lb = [];
params0 = [];
sampParam0 = @(lb, ub, n) rand(1,n) .* (ub-lb) + lb;

nparams_full = length(paramInd);
nparamsAll = sum(nLocTrue.^paramInd);
% organize ub, lb and param0 given a certain param combination
rng default % For reproducibility
for iparam = 1:nparams_full
    if paramInd(iparam) == 1
        ub = [ub, ones(1, nLocTrue)*ub_full(iparam)];
        lb = [lb, ones(1,nLocTrue)*lb_full(iparam)];
        params0 = [params0, sampParam0(lb_full(iparam), ub_full(iparam), nLocTrue)]; % set param0
    else
        ub = [ub, ub_full(iparam)];
        lb = [lb, lb_full(iparam)];
        params0 = [params0, sampParam0(lb_full(iparam), ub_full(iparam), 1)]; % set param0
    end
end

assert(length(ub) == nparamsAll)
assert(length(lb) == nparamsAll)
assert(length(params0) == nparamsAll)
