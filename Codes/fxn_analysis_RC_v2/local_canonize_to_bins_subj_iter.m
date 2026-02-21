function Xout = local_canonize_to_bins_subj_iter(Xin, nBins, nSubj)
% Force Xin into [nBins x nSubj x nIter] by identifying the bin/subj dims.
X = Xin;
sz = size(X);
nd = ndims(X);

% pad to 3D if needed
if nd == 2
    % could be [nBins x nSubj] (no iter)
    X = reshape(X, sz(1), sz(2), 1);
    sz = size(X);
    nd = 3;
elseif nd > 3
    % collapse trailing dims into iter
    X = reshape(X, sz(1), sz(2), prod(sz(3:end)));
    sz = size(X);
    nd = 3;
end

% Identify dims
dims = 1:nd;
iBin  = find(sz == nBins, 1, 'first');
iSubj = find(sz == nSubj, 1, 'first');

if isempty(iBin) || isempty(iSubj) || iBin == iSubj
    error('Cannot canonize: expected one dim==nBins and one dim==nSubj. Got size %s.', mat2str(sz));
end

iIter = setdiff(dims, [iBin iSubj]);
if isempty(iIter), iIter = 3; end

Xout = permute(X, [iBin iSubj iIter]);
end