noiseCST = .2;
gaborCST = .3;
nTrials = 4e3;
nIter = 10;
cSDT_true = 0;
lambda_whiten = 0;

iModelB_sim = 6;
Nmul_true = .5;
for Nadd_true = [5]
    Nshared_true = 1;
    % OOD_sim(noiseCST, gaborCST, nTrials, Nmul_true, Nadd_true, Nshared_true, cSDT_true, lambda_whiten, iModelB_sim, nIter)
end
