
gaborCST=.5;
nTrials=1e4;
nIter = 20;
cSDT_true = 0;

noiseCST=.2;
lambda_whiten = 0;

for Nmul_true = [.6 .9]
    for Nadd_true = [6,9]
        for Nshared_true = [6,9]
            for iModelB_sim = 1:7
                OOD_sim(noiseCST, gaborCST, nTrials, Nmul_true, Nadd_true, Nshared_true, cSDT_true, ...
                    lambda_whiten, iModelB_sim, nIter)
            end
        end
    end
end
