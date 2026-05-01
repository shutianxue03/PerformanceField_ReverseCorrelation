
gaborCST=.5;
nTrials=1e4;
nIter = 20;
cSDT_true = 0;

noiseCST=.2;
flag_regressType = 2;
flag_incluCrit = 1;
C_contribution = 0;
lambda_whiten = 0;
nBasisORI = 6;
nBasisSF = 5;

for Nmul_true = .6
    for Nadd_true = 6
        for Nshared_true = 6
            for iModelB_sim = 1:7
                OOD_sim(noiseCST, gaborCST, nTrials, Nmul_true, Nadd_true, Nshared_true, cSDT_true, ...
                    lambda_whiten, flag_regressType, flag_incluCrit, C_contribution, iModelB_sim, nIter, nBasisORI, nBasisSF)
            end
        end
    end
end
