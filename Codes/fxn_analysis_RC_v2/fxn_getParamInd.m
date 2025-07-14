function indParamIncl = fxn_getParamInd(nParams_full)
nCand = 2^nParams_full;
indParamIncl = zeros(nCand, nParams_full);
im_true = 1;
nmodels_all = size(1,nParams_full);

for iParam = 1:nParams_full
    comb = nchoosek(1:nParams_full, iParam);
    nmodels_all(iParam) = size(comb,1);
    nchosen = size(comb,2);
    for im = 1:nmodels_all(iParam)
        for ic = 1:nchosen, indParamIncl(im_true, comb(im, ic)) = 1; end
        im_true = im_true+1;
    end
end

assert(sum(nmodels_all) == nCand-1)
indParamIncl = (indParamIncl==1);
end