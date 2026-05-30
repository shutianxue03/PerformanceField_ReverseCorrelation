
noiseCST = .2;
gaborCST = .3;
nTrials = 4e3;
nIter = 24;
cSDT_true = 0;


for iModelB_sim = [5] % 5=Nmul only; 6=Nadd only; 7=Nshared only

    for Nmul_true = .5
        for Nadd_true = [1]
            for Nshared_true = 0
                OOD_sim(noiseCST, gaborCST, nTrials, Nmul_true, Nadd_true, Nshared_true, cSDT_true, iModelB_sim, nIter)
            end
        end
    end
end