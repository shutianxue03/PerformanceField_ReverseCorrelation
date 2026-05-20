noiseCST = .2;
gaborCST = .3;
nTrials = 1e4;
nIter = 20;
cSDT_true = 0;
lambda_whiten = 0;
nBasisORI = 6;
basisWidthORI = 0.9;
nBasisSF = 6;
basisWidthSF = 0.6;
Nmul_true=.5;
Nadd_true=1;
Nshared_true=1;
iModelB_sim = 6;

OOD_sim(noiseCST, gaborCST, nTrials, Nmul_true, Nadd_true, Nshared_true, cSDT_true, lambda_whiten, iModelB_sim, nIter, nBasisORI, basisWidthORI, nBasisSF, basisWidthSF)


