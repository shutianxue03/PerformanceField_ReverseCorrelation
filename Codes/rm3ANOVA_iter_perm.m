function stats = rm3ANOVA_iter_perm(nLL_allCond, nPerm, CI_level)
% fastRM3way_ANOVA_iter_perm
% 3-way repeated-measures ANOVA (all factors within-subject) for
% nLL_allCond: [nA x nB x nL x nSubj x nIter]
%
% Outputs:
%   stats.F_allIter:      [nIter x 7]  (A,B,L,AB,AL,BL,ABL)
%   stats.eta2_allIter:   [nIter x 7]
%   stats.F_med/lb/ub:    [1 x 7]
%   stats.eta2_med/lb/ub: [1 x 7]
%   stats.p_perm:         [1 x 7] permutation p-values (right-tailed for F)

if nargin < 2 || isempty(nPerm), nPerm = 1e4; end
if nargin < 3 || isempty(CI_level), CI_level = 0.95; end

% ---------- validate ----------
assert(ndims(nLL_allCond) == 5, 'nLL_allCond must be 5D: [A x B x L x Subj x Iter].');
[nA,nB,nL,nS,nIter] = size(nLL_allCond);

% ---------- per-iteration stats ----------
F_allIter    = nan(nIter, 7);
eta2_allIter = nan(nIter, 7);

parfor iIter = 1:nIter
    D = nLL_allCond(:,:,:,:,iIter);               % [A B L S]
    [F, eta2] = rm3way_ss_fast(D);                % both 1x7
    F_allIter(iIter,:)    = F;
    eta2_allIter(iIter,:) = eta2;
end

% median + CI across iterations
[F_med,   F_lb,   F_ub]   = getCI(F_allIter,    1, 1, CI_level);
[eta_med, eta_lb, eta_ub] = getCI(eta2_allIter, 1, 1, CI_level);

% ---------- permutation on subject-level medians across iterations ----------
% observer-level point estimate = median across iterations
Dmed = getCI(nLL_allCond, 1, ndims(nLL_allCond));   % median over Iter -> [A B L S]
% NOTE: your getCI returns [nS x ...] sometimes; enforce shape
Dmed = reshape(Dmed, [nA,nB,nL,nS]);

[F_obs, ~] = rm3way_ss_fast(Dmed);  % 1x7

% Pre-vectorize condition dimension for fast within-subject shuffles
nCond = nA*nB*nL;
Dmed_vec = reshape(Dmed, [nCond, nS]);  % [cond x subj]

F_perm = nan(nPerm, 7);

% Use parfor if you want (usually faster for nPerm>=5000)
parfor iPerm = 1:nPerm
    X = Dmed_vec;

    % Shuffle condition labels WITHIN subject:
    % null = exchangeability across all condition cells within each subject
    for s = 1:nS
        X(:,s) = X(randperm(nCond), s);
    end

    Dp = reshape(X, [nA,nB,nL,nS]);
    Fp = rm3way_ss_fast(Dp);         % 1x7
    F_perm(iPerm,:) = Fp;
end

% Right-tailed p-values for F (F is nonnegative; “two-tailed F” isn’t meaningful)
p_perm = (1 + sum(F_perm >= F_obs, 1)) ./ (nPerm + 1);

% ---------- pack ----------
stats = struct();
stats.effectNames = {'A','B','L','AB','AL','BL','ABL'};
stats.F_med = F_med;   
stats.F_lb = F_lb;   
stats.F_ub = F_ub;
stats.eta2_med = eta_med; 
stats.eta2_lb = eta_lb; 
stats.eta2_ub = eta_ub;
stats.p_perm = round(p_perm, 3);

end


function [F, eta2p] = rm3way_ss_fast(D)
% rm3way_ss_fast
% Fast 3-way RM ANOVA SS-based computation for balanced design.
% D: [A B L S] numeric, may include NaNs (handled via 'omitnan' means).

[nA,nB,nL,nS] = size(D);

% Means
G   = mean(D, 'all', 'omitnan');                 % grand mean
Ms  = mean(D, [1 2 3], 'omitnan');               % [1 1 1 S]
Ma  = mean(D, [2 3 4], 'omitnan');               % [A 1 1 1]
Mb  = mean(D, [1 3 4], 'omitnan');               % [1 B 1 1]
Ml  = mean(D, [1 2 4], 'omitnan');               % [1 1 L 1]

Mab = mean(D, [3 4], 'omitnan');                 % [A B 1 1]
Mal = mean(D, [2 4], 'omitnan');                 % [A 1 L 1]
Mbl = mean(D, [1 4], 'omitnan');                 % [1 B L 1]
Mabl= mean(D, 4, 'omitnan');                     % [A B L 1]

% Subject x factor marginal means (for error terms)
Mas  = mean(D, [2 3], 'omitnan');                % [A 1 1 S]
Mbs  = mean(D, [1 3], 'omitnan');                % [1 B 1 S]
Mls  = mean(D, [1 2], 'omitnan');                % [1 1 L S]

Mabs = mean(D, 3, 'omitnan');                    % [A B 1 S]
Mals = mean(D, 2, 'omitnan');                    % [A 1 L S]
Mbls = mean(D, 1, 'omitnan');                    % [1 B L S]

% --- SS for effects (within-subject) ---
SSa   = nS*nB*nL * sum((Ma - G).^2, 'all', 'omitnan');
SSb   = nS*nA*nL * sum((Mb - G).^2, 'all', 'omitnan');
SSl   = nS*nA*nB * sum((Ml - G).^2, 'all', 'omitnan');

SSab  = nS*nL * sum((Mab - Ma - Mb + G).^2, 'all', 'omitnan');
SSal  = nS*nB * sum((Mal - Ma - Ml + G).^2, 'all', 'omitnan');
SSbl  = nS*nA * sum((Mbl - Mb - Ml + G).^2, 'all', 'omitnan');

SSabl = nS * sum((Mabl ...
    - Mab - Mal - Mbl ...
    + Ma + Mb + Ml - G).^2, 'all', 'omitnan');

% --- SS error terms = subject x effect interactions ---
SSas   = nB*nL * sum((Mas - Ma - Ms + G).^2, 'all', 'omitnan');
SSbs   = nA*nL * sum((Mbs - Mb - Ms + G).^2, 'all', 'omitnan');
SSls   = nA*nB * sum((Mls - Ml - Ms + G).^2, 'all', 'omitnan');

SSabs  = nL * sum((Mabs ...
    - Mab - Mas - Mbs ...
    + Ma + Mb + Ms - G).^2, 'all', 'omitnan');

SSals  = nB * sum((Mals ...
    - Mal - Mas - Mls ...
    + Ma + Ml + Ms - G).^2, 'all', 'omitnan');

SSbls  = nA * sum((Mbls ...
    - Mbl - Mbs - Mls ...
    + Mb + Ml + Ms - G).^2, 'all', 'omitnan');

SSabls = sum((D ...
    - Mabl - Mabs - Mals - Mbls ...
    + Mab + Mal + Mbl + Mas + Mbs + Mls ...
    - Ma - Mb - Ml - Ms + G).^2, 'all', 'omitnan');

% dfs
dfa   = nA - 1;
dfb   = nB - 1;
dfl   = nL - 1;

dfas  = (nS - 1) * dfa;
dfbs  = (nS - 1) * dfb;
dfls  = (nS - 1) * dfl;

dfab  = dfa * dfb;
dfal  = dfa * dfl;
dfbl  = dfb * dfl;

dfabs = (nS - 1) * dfab;
dfals = (nS - 1) * dfal;
dfbls = (nS - 1) * dfbl;

dfabl  = dfa * dfb * dfl;
dfabls = (nS - 1) * dfabl;

% MS and F
MSa   = SSa   / dfa;    MSas   = SSas   / dfas;    Fa   = MSa   / MSas;
MSb   = SSb   / dfb;    MSbs   = SSbs   / dfbs;    Fb   = MSb   / MSbs;
MSl   = SSl   / dfl;    MSls   = SSls   / dfls;    Fl   = MSl   / MSls;

MSab  = SSab  / dfab;   MSabs  = SSabs  / dfabs;   Fab  = MSab  / MSabs;
MSal  = SSal  / dfal;   MSals  = SSals  / dfals;   Fal  = MSal  / MSals;
MSbl  = SSbl  / dfbl;   MSbls  = SSbls  / dfbls;   Fbl  = MSbl  / MSbls;

MSabl = SSabl / dfabl;  MSabls = SSabls / dfabls;  Fabl = MSabl / MSabls;

F = [Fa Fb Fl Fab Fal Fbl Fabl];

% partial eta^2 for RM: SS_effect / (SS_effect + SS_SxEffect)
eta2p = [ ...
    SSa  /(SSa  + SSas), ...
    SSb  /(SSb  + SSbs), ...
    SSl  /(SSl  + SSls), ...
    SSab /(SSab + SSabs), ...
    SSal /(SSal + SSals), ...
    SSbl /(SSbl + SSbls), ...
    SSabl/(SSabl+ SSabls) ...
    ];
end