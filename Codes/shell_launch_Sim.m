noiseCST = .2;

gaborCST=.3;
nTrials=1e4;
nIter = 16;
cSDT_true = 0;
lambda_whiten = 0;
iModelB_sim = 6;

for flag_whitenDV = 1%0:1
    switch flag_whitenDV
        case 0, Nmul_true = .5; Nadd_true = 5; Nshared_true = 5;
            case 1, Nmul_true = .2; Nadd_true = .2; Nshared_true = .2;
    end

    OOD_sim(noiseCST, gaborCST, nTrials, Nmul_true, Nadd_true, Nshared_true, cSDT_true, lambda_whiten, flag_whitenDV, iModelB_sim, nIter)

end
