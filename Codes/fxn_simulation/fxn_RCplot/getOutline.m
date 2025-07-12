function outline = getOutline(nORI, nSF, data2D, nB, pCriterion)
% get outline
if ndims(data2D) == 2, error('ALERT:The input matrix is not bootstrapped'), end
if size(data2D,2) ~= nORI, error('ALERT:The shape of input matrix is wrong'), end
p_all = nan(nORI, nSF);
for iORI = 1:nORI
    for iSF = 1:nSF
        p_all(iORI, iSF) = sum(squeeze(data2D(:, iORI, iSF)) > 0)/nB;
    end
end
outline = p_all >= 1-pCriterion;
end
