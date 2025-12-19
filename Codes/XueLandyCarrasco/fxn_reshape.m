function A_out = fxn_reshape(A_in, nSubj, nCond, varName)
% Identify subject/condition dims by size matching, then permute to:
%   [nBoot x nSubj x nCond]
% Failure mode (rare): if nSubj == nCond, size-based matching is ambiguous.
sz = size(A_in);
indSubj = find(sz == nSubj, 1, 'first');
indCond = find(sz == nCond, 1, 'first');

assert(~isempty(indSubj), '%s: cannot find a dimension matching nSubj=%d. size=%s', varName, nSubj, mat2str(sz));
assert(~isempty(indCond), '%s: cannot find a dimension matching nCond=%d. size=%s', varName, nCond, mat2str(sz));
assert(indSubj ~= indCond, '%s: ambiguous dims (nSubj and nCond match same dim). size=%s', varName, mat2str(sz));

indBoot = setdiff(1:3, [indSubj, indCond]);
assert(numel(indBoot)==1, '%s: cannot identify boot dim uniquely. size=%s', varName, mat2str(sz));

A_out = permute(A_in, [indBoot, indSubj, indCond]);
end