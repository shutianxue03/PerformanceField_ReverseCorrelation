noiseCST = .2;
gaborCST = .3;
nTrials = 4e3;
nIter = 20;
cSDT_true = 0;
lambda_whiten = 0;

nBasisORI = 3;
basisWidthORI = 0.9;

for nBasisSF = 1:2:9
for basisWidthSF = .1:.2:.9

iModelB_sim = 6;
Nmul_true = .5;
Nadd_true = 5;
Nshared_true = 1;

OOD_sim(noiseCST, gaborCST, nTrials, Nmul_true, Nadd_true, Nshared_true, cSDT_true, lambda_whiten, iModelB_sim, nIter, nBasisORI, basisWidthORI, nBasisSF, basisWidthSF)
end
end
