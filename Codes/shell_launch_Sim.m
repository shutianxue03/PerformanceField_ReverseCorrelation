
noiseCST = .2;
gaborCST = .3;
nTrials = 4e3;
nIter = 10;
cSDT_true = 0;

iModelB_sim = 6;

Nmul_true = .5;
for Nadd_true = [1]
    Nshared_true = 1;
    OOD_sim(noiseCST, gaborCST, nTrials, Nmul_true, Nadd_true, Nshared_true, cSDT_true, iModelB_sim, nIter)
end
