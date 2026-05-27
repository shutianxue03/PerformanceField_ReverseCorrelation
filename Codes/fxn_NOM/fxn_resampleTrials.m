function R = fxn_resampleTrials(dataMatrix, e3D_allT, iSess_select, iLoc_all, ratio_split)
% Resample trials by pair ID, then split pairs into full/template/train/test sets.
%
% Inputs
%   dataMatrix   : trial metadata matrix.
%                  Uses columns: 2=session, 5=location, 6=iPRS, 8=iPair,
%                  9=resp, 10=RT, 11=contrast.
%   e3D_allT     : trial-by-orientation-by-SF energy tensor.
%   iSess_select : sessions to include.
%   iLoc_all     : location index (scalar for this function version).
%   ratio_split  : [template train test] fractions summing to 1.
%
% Output
%   R            : struct with fields for each split
%                  (e3D_*, iPRS_*, resp_*, cst_*, respC_*, iPair_*, RT_*).

% FULL set (no split) - optional
e3D_full_rand = [];
iPRS_full_rand = [];
resp_full_rand = [];
cst_full_rand = [];
respC_full_rand = [];
iPair_full_rand = [];
RT_full_rand = [];

% TEMPLATE set (for RC/template derivation)
e3D_tmpl_rand = [];
iPRS_tmpl_rand = [];
resp_tmpl_rand = [];
cst_tmpl_rand = [];
respC_tmpl_rand = [];
iPair_tmpl_rand = [];
RT_tmpl_rand = [];

% Training set (for parameter estimation)
e3D_train_rand = [];
iPRS_train_rand = [];
resp_train_rand = [];
cst_train_rand = [];
respC_train_rand = [];
iPair_train_rand = [];
RT_train_rand = [];

% Test set (for log likelihood + predictions)
e3D_test_rand = [];
iPRS_test_rand = [];
resp_test_rand = [];
cst_test_rand = [];
respC_test_rand = [];
iPair_test_rand = [];
RT_test_rand = [];

iLoc = iLoc_all;

% Collect data for this location
e3D_nonrand = [];
iPRS_nonrand = [];
iPair_nonrand = [];
resp_nonrand = [];
cst_nonrand = [];
RT_nonrand = [];

for iSS = 1:length(iSess_select)
    indLoc = (dataMatrix(:, 5) == iLoc) & (dataMatrix(:, 2) == iSess_select(iSS));

    e3D_nonrand = cat(1, e3D_nonrand, e3D_allT(indLoc, :, :));
    iPRS_nonrand = [iPRS_nonrand; dataMatrix(indLoc, 6)];
    iPair_nonrand = [iPair_nonrand; dataMatrix(indLoc, 8)];
    resp_nonrand = [resp_nonrand; dataMatrix(indLoc, 9)];
    cst_nonrand = [cst_nonrand; dataMatrix(indLoc, 11)];
    RT_nonrand = [RT_nonrand; dataMatrix(indLoc, 10)];
end

% Unique pair IDs define sampling units (prevents splitting a pair across sets).
iPair_unik_nonrand = unique(iPair_nonrand);
nPairs_total = numel(iPair_unik_nonrand);
assert(nPairs_total > 0, 'No pairs found for this location/session selection.');

% Randomize pairs, then split pairs into template/train/test
indUnikPair_rand = randperm(nPairs_total, nPairs_total);

nPairs_tmpl = floor(ratio_split(1) * nPairs_total);
nPairs_train = floor(ratio_split(2) * nPairs_total);
nPairs_test = nPairs_total - nPairs_tmpl - nPairs_train;

% Guardrails: keep all three splits non-empty when possible.
nPairs_tmpl = max(nPairs_tmpl, 1);
nPairs_train = max(nPairs_train, 1);
nPairs_test = nPairs_total - nPairs_tmpl - nPairs_train;
if nPairs_test < 1
    [~, idxMax] = max([nPairs_tmpl, nPairs_train]);
    if idxMax == 1
        nPairs_tmpl = nPairs_tmpl - 1;
    else
        nPairs_train = nPairs_train - 1;
    end
    nPairs_test = 1;
end
assert(nPairs_tmpl + nPairs_train + nPairs_test == nPairs_total)

indUnikPair_rand_tmpl_ = indUnikPair_rand(1:nPairs_tmpl);
indUnikPair_rand_train_ = indUnikPair_rand(nPairs_tmpl+1 : nPairs_tmpl+nPairs_train);
indUnikPair_rand_test_ = indUnikPair_rand(nPairs_tmpl+nPairs_train+1 : end);

% Build FULL set from all randomized pairs.
[e3D_full_, iPRS_full_, resp_full_, cst_full_, respC_full_, iPair_full_, RT_full_] = fxn_getRespC( ...
    e3D_nonrand, indUnikPair_rand, iPair_nonrand, iPair_unik_nonrand, ...
    iPRS_nonrand, resp_nonrand, cst_nonrand, RT_nonrand);

e3D_full_rand = cat(1, e3D_full_rand, e3D_full_);
iPRS_full_rand = [iPRS_full_rand; iPRS_full_];
resp_full_rand = [resp_full_rand; resp_full_];
cst_full_rand = [cst_full_rand; cst_full_];
respC_full_rand = [respC_full_rand; respC_full_];
iPair_full_rand = [iPair_full_rand; iPair_full_];
RT_full_rand = [RT_full_rand; RT_full_];

% Build TEMPLATE set from template split pairs.
[e3D_tmpl_, iPRS_tmpl_, resp_tmpl_, cst_tmpl_, respC_tmpl_, iPair_tmpl_, RT_tmpl_] = fxn_getRespC( ...
    e3D_nonrand, indUnikPair_rand_tmpl_, iPair_nonrand, iPair_unik_nonrand, ...
    iPRS_nonrand, resp_nonrand, cst_nonrand, RT_nonrand);

e3D_tmpl_rand = cat(1, e3D_tmpl_rand, e3D_tmpl_);
iPRS_tmpl_rand = [iPRS_tmpl_rand; iPRS_tmpl_];
resp_tmpl_rand = [resp_tmpl_rand; resp_tmpl_];
cst_tmpl_rand = [cst_tmpl_rand; cst_tmpl_];
respC_tmpl_rand = [respC_tmpl_rand; respC_tmpl_];
iPair_tmpl_rand = [iPair_tmpl_rand; iPair_tmpl_];
RT_tmpl_rand = [RT_tmpl_rand; RT_tmpl_];

% Build TRAIN set from train split pairs.
[e3D_train_, iPRS_train_, resp_train_, cst_train_, respC_train_, iPair_train_, RT_train_] = fxn_getRespC( ...
    e3D_nonrand, indUnikPair_rand_train_, iPair_nonrand, iPair_unik_nonrand, ...
    iPRS_nonrand, resp_nonrand, cst_nonrand, RT_nonrand);

e3D_train_rand = cat(1, e3D_train_rand, e3D_train_);
iPRS_train_rand = [iPRS_train_rand; iPRS_train_];
resp_train_rand = [resp_train_rand; resp_train_];
cst_train_rand = [cst_train_rand; cst_train_];
respC_train_rand = [respC_train_rand; respC_train_];
iPair_train_rand = [iPair_train_rand; iPair_train_];
RT_train_rand = [RT_train_rand; RT_train_];

% Build TEST set from test split pairs.
[e3D_test_, iPRS_test_, resp_test_, cst_test_, respC_test_, iPair_test_, RT_test_] = fxn_getRespC( ...
    e3D_nonrand, indUnikPair_rand_test_, iPair_nonrand, iPair_unik_nonrand, ...
    iPRS_nonrand, resp_nonrand, cst_nonrand, RT_nonrand);

e3D_test_rand = cat(1, e3D_test_rand, e3D_test_);
iPRS_test_rand = [iPRS_test_rand; iPRS_test_];
resp_test_rand = [resp_test_rand; resp_test_];
cst_test_rand = [cst_test_rand; cst_test_];
respC_test_rand = [respC_test_rand; respC_test_];
iPair_test_rand = [iPair_test_rand; iPair_test_];
RT_test_rand = [RT_test_rand; RT_test_];

% Package all outputs into one struct for parfor-safe passing.
R = struct();
R.e3D_full_rand = e3D_full_rand;
R.iPRS_full_rand = iPRS_full_rand;
R.resp_full_rand = resp_full_rand;
R.cst_full_rand = cst_full_rand;
R.respC_full_rand = respC_full_rand;
R.iPair_full_rand = iPair_full_rand;
R.RT_full_rand = RT_full_rand;
R.e3D_tmpl_rand = e3D_tmpl_rand;
R.iPRS_tmpl_rand = iPRS_tmpl_rand;
R.resp_tmpl_rand = resp_tmpl_rand;
R.cst_tmpl_rand = cst_tmpl_rand;
R.respC_tmpl_rand = respC_tmpl_rand;
R.iPair_tmpl_rand = iPair_tmpl_rand;
R.RT_tmpl_rand = RT_tmpl_rand;
R.e3D_train_rand = e3D_train_rand;
R.iPRS_train_rand = iPRS_train_rand;
R.resp_train_rand = resp_train_rand;
R.cst_train_rand = cst_train_rand;
R.respC_train_rand = respC_train_rand;
R.iPair_train_rand = iPair_train_rand;
R.RT_train_rand = RT_train_rand;
R.e3D_test_rand = e3D_test_rand;
R.iPRS_test_rand = iPRS_test_rand;
R.resp_test_rand = resp_test_rand;
R.cst_test_rand = cst_test_rand;
R.respC_test_rand = respC_test_rand;
R.iPair_test_rand = iPair_test_rand;
R.RT_test_rand = RT_test_rand;

end