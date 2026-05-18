noiseCST = .2;
gaborCST = .3;
nTrials = 1e4;
nIter = 100;
cSDT_true = 0;
lambda_whiten = 0;

for iModelB_sim = 6
    Nmul_true = .5; Nadd_true = 5; Nshared_true = 5;
    OOD_sim(noiseCST, gaborCST, nTrials, Nmul_true, Nadd_true, Nshared_true, cSDT_true, lambda_whiten, iModelB_sim, nIter)
end