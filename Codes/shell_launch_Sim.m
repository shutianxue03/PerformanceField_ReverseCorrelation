noiseCST=.2;
gaborCST=.2;
nTrials=1e4;
Nmul_true=0;
Nadd_true=0;
Nshared_true=0;
flag_incluCrit = 1;
iModelA_sim = 1;
iModelB_sim = 1;
nIter = 10;
flag_regressType = 2;
lambda_whiten = .5;
C_contribution = 0;
iModelA_sim=1

Nadd_true = 0;
for Nmul_true = [.2]
    for iModelB_sim = [4,1]
        for gaborCST=[.3];
            for Cz_true = 0%[-.4:.2:.2]
                for lambda_whiten = 1

                    % OOD_sim(.2, .2, 4e3, 0, 0, 0, Cz_true, lambda_whiten, flag_regressType, iModelA_sim, 1, 10)
                    OOD_sim(noiseCST, gaborCST, nTrials, Nmul_true, Nadd_true, Nshared_true, Cz_true, ...
                        lambda_whiten, flag_regressType, flag_incluCrit, C_contribution, iModelA_sim, iModelB_sim, nIter)
                end
            end
        end
    end
end

Nmul_true = 0;
for Nadd_true = [5]
    for iModelB_sim = [3,1]
        for gaborCST=[.3];
            for Cz_true = 0%[-.4:.2:.2]
                for lambda_whiten = 1
                    % OOD_sim(.2, .2, 4e3, 0, 0, 0, Cz_true, lambda_whiten, flag_regressType, iModelA_sim, 1, 10)
                    OOD_sim(noiseCST, gaborCST, nTrials, Nmul_true, Nadd_true, Nshared_true, Cz_true, ...
                        lambda_whiten, flag_regressType, flag_incluCrit, C_contribution, iModelA_sim, iModelB_sim, nIter)
                end
            end
        end
    end
end

