
flag_stage = 2;
nBoot = 20;

for isubj = 15:-1:1
    for iLocComb = 1:5
        for iModelA = 1:2
            switch flag_stage
                case 1
                    OOD_NOM_Trialwise_compIV(isubj, iLocComb, iModelA, nBoot);
                case 2
                    for iModelB = 1:2
                        OOD_NOM_Trialwise_fitNOM(isubj, iLocComb, iModelA, iModelB, nBoot);
                    end
            end
        end
    end
end
