noiseCST=.2;
gaborCST=.4;
nTrials=1e4;
flag_incluCrit = 1;
nIter = 2;
flag_regressType = 2;
C_contribution = 0;
lambda_whiten = 1;
Cz_true = 0;

for Nmul_true = .1
    for Nadd_true = 10
        for Nshared_true = 10
            for iModelB_sim = 1:4
                OOD_sim(noiseCST, gaborCST, nTrials, Nmul_true, Nadd_true, Nshared_true, Cz_true, ...
                    lambda_whiten, flag_regressType, flag_incluCrit, C_contribution, iModelB_sim, nIter)
            end
        end
    end
end
