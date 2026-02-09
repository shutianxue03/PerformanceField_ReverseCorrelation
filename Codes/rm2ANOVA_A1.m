function [stats, str_stats] = rm2ANOVA_A1(nLL_allCond, nPerm, CI_level)
% Two-way repeated-measures ANOVA for ModelB x Loc, computed per iteration.
%
% INPUT
%   nLL_allCond : [nB x nL x nSubj x nIter]
%   nPerm       : e.g., 1e4
%   CI_level    : e.g., 0.95
%
% OUTPUT (struct)
%   stats.F_allIter, stats.eta2p_allIter : per-iteration stats
%   stats.F_med/lb/ub, stats.eta2p_med/lb/ub : across-iteration summaries
%   stats.p_perm : permutation p-values (right-tailed) for each effect

if nargin < 2 || isempty(CI_level), CI_level = 0.95; end
if nargin < 3 || isempty(nPerm),    nPerm    = 1e4;  end

assert(ndims(nLL_allCond) == 4, 'Input must be [nB x nL x nSubj x nIter].');

[nB, nL, nS, nIter] = size(nLL_allCond);

% --- prealloc ---
F_B_allIter  = nan(nIter,1);
F_L_allIter  = nan(nIter,1);
F_BL_allIter = nan(nIter,1);

eta2p_B_allIter  = nan(nIter,1);
eta2p_L_allIter  = nan(nIter,1);
eta2p_BL_allIter = nan(nIter,1);

% =========================================================
% (1) Per-iteration RM two-way ANOVA: fast SS partition
% =========================================================
for iIter = 1:nIter
    X = nLL_allCond(:,:,:,iIter);  % [nB x nL x nS]
    if all(isnan(X(:))), continue; end

    [F, eta2p] = rm_twoway_fast(X);

    F_B_allIter(iIter)     = F.B;
    F_L_allIter(iIter)     = F.Loc;
    F_BL_allIter(iIter)    = F.BxLoc;

    eta2p_B_allIter(iIter)  = eta2p.B;
    eta2p_L_allIter(iIter)  = eta2p.Loc;
    eta2p_BL_allIter(iIter) = eta2p.BxLoc;
end

% --- CI across iterations (median + CI) ---
[F_B_med,  F_B_lb,  F_B_ub]  = getCI(F_B_allIter,  1,1, CI_level);
[F_L_med,  F_L_lb,  F_L_ub]  = getCI(F_L_allIter,  1,1, CI_level);
[F_BL_med, F_BL_lb, F_BL_ub] = getCI(F_BL_allIter, 1,1, CI_level);

[eta_B_med,  eta_B_lb,  eta_B_ub]  = getCI(eta2p_B_allIter,  1,1, CI_level);
[eta_L_med,  eta_L_lb,  eta_L_ub]  = getCI(eta2p_L_allIter,  1,1, CI_level);
[eta_BL_med, eta_BL_lb, eta_BL_ub] = getCI(eta2p_BL_allIter, 1,1, CI_level);

% =========================================================
% (2) Permutation test on point estimates (median across iters)
%     Permute condition labels WITHIN subject across the nB*nL cells.
% =========================================================
X_med = median(nLL_allCond, 4, 'omitnan');  % [nB x nL x nS]

[F_obs, ~] = rm_twoway_fast(X_med);

Fperm_B  = nan(nPerm,1);
Fperm_L  = nan(nPerm,1);
Fperm_BL = nan(nPerm,1);

nCond = nB*nL;

parfor iPerm = 1:nPerm
    Xp = X_med; % [nB x nL x nS]

    % permute cell labels within each subject
    for s = 1:nS
        v = reshape(Xp(:,:,s), [nCond 1]);
        if all(isnan(v)), continue; end
        rp = randperm(nCond);
        v = v(rp);
        Xp(:,:,s) = reshape(v, [nB nL]);
    end

    Fp = rm_twoway_fast(Xp);

    Fperm_B(iPerm)  = Fp.B;
    Fperm_L(iPerm)  = Fp.Loc;
    Fperm_BL(iPerm) = Fp.BxLoc;
end

% right-tailed permutation p-values (F >= 0)
p_B  = (1 + sum(Fperm_B  >= F_obs.B))     / (nPerm + 1);
p_L  = (1 + sum(Fperm_L  >= F_obs.Loc))   / (nPerm + 1);
p_BL = (1 + sum(Fperm_BL >= F_obs.BxLoc)) / (nPerm + 1);

% =========================================================
% Pack outputs
% =========================================================

stats = struct();
stats.effectNames = {'B','L','BL'};
stats.F_med = [F_B_med, F_L_med, F_BL_med];   
stats.F_lb = [F_B_lb, F_L_lb, F_BL_lb];   
stats.F_ub = [F_B_ub, F_L_ub, F_BL_ub];
stats.eta2_med = [eta_B_med, eta_L_med, eta_BL_med]; 
stats.eta2_lb = [eta_B_lb, eta_L_lb, eta_BL_lb]; 
stats.eta2_ub = [eta_B_ub, eta_L_ub, eta_BL_ub];
stats.p_perm = round([p_B, p_L, p_BL], 3);

% =========================================================
% Create string
% =========================================================
nEffects = numel(stats.effectNames);
str_stats = [];

for iEffect = 1:nEffects
    str_stats = [str_stats, sprintf( ...
        '%s: F=%.2f, p=%.3f, eta2=%.2f [%.2f, %.2f]\n', ...
        stats.effectNames{iEffect}, ...
        stats.F_med(iEffect), ...
        stats.p_perm(iEffect), ...
        stats.eta2_med(iEffect), ...
        stats.eta2_lb(iEffect), ...
        stats.eta2_ub(iEffect))];
end
end

%%
function [F, eta2p] = rm_twoway_fast(X)
% Fast two-way repeated-measures ANOVA for within-subject factors:
%   Factor 1: B (nB levels)
%   Factor 2: Loc (nL levels)
%
% INPUT: X is [nB x nL x nS]
% OUTPUT:
%   F fields: .B, .Loc, .BxLoc
%   eta2p fields: partial eta^2 for each effect

[nB, nL, nS] = size(X);

% means
GM  = mean(X(:), 'omitnan');
mS  = squeeze(mean(mean(X,1,'omitnan'),2,'omitnan'));      % [nS x 1]
mB  = squeeze(mean(mean(X,2,'omitnan'),3,'omitnan'));      % [nB x 1]
mL  = squeeze(mean(mean(X,1,'omitnan'),3,'omitnan'));      % [nL x 1]
mBL = squeeze(mean(X,3,'omitnan'));                        % [nB x nL]
mBS = squeeze(mean(X,2,'omitnan'));                        % [nB x nS]
mLS = squeeze(mean(X,1,'omitnan'));                        % [nL x nS]

% Force consistent shapes
mB = mB(:);          % nB x 1
mL = mL(:);          % nL x 1
mS = mS(:);          % nS x 1

% SS partitions (balanced RM; assumes missingness is limited—ok for your pipeline)
SS_T  = sum((X(:) - GM).^2);
SS_S  = nB*nL * sum((mS - GM).^2);
SS_B  = nS*nL * sum((mB - GM).^2);
SS_L  = nS*nB * sum((mL - GM).^2);

% 2-way interaction SS terms
SS_BL = nS * sum( (mBL - mB - mL.' + GM).^2 , 'all');   % (nB x nL)
SS_BS = nL * sum( (mBS - mB - mS.' + GM).^2 , 'all');   % (nB x nS)
SS_LS = nB * sum( (mLS - mL - mS.' + GM).^2 , 'all');   % (nL x nS)

SS_BLS = SS_T - SS_S - SS_B - SS_L - SS_BL - SS_BS - SS_LS;

% dfs
df_S   = nS - 1;
df_B   = nB - 1;
df_L   = nL - 1;
df_BL  = (nB - 1)*(nL - 1);

df_BS  = df_B * df_S;
df_LS  = df_L * df_S;
df_BLS = df_BL * df_S;

% MS + F
MS_B   = SS_B   / df_B;
MS_L   = SS_L   / df_L;
MS_BL  = SS_BL  / df_BL;

MS_BS  = SS_BS  / df_BS;
MS_LS  = SS_LS  / df_LS;
MS_BLS = SS_BLS / df_BLS;

F = struct();
F.B     = MS_B  / MS_BS;
F.Loc   = MS_L  / MS_LS;
F.BxLoc = MS_BL / MS_BLS;

% partial eta^2 using the effect’s error term
eta2p = struct();
eta2p.B     = SS_B  / (SS_B  + SS_BS);
eta2p.Loc   = SS_L  / (SS_L  + SS_LS);
eta2p.BxLoc = SS_BL / (SS_BL + SS_BLS);

end