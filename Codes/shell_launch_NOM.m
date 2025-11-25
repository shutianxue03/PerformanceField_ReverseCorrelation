
flag_stage = 2;
nBoot = 10;

for isubj = [1,3]
    for iLocComb = [1, 8, 6, 7]
        for iModelA = 1
            switch flag_stage
                case 1
                    OOD_NOM_Trialwise_compIV(isubj, iLocComb, iModelA, nBoot);
                case 2
                    for iModelB = [1,2]
                        OOD_NOM_Trialwise_fitNOM(isubj, iLocComb, iModelA, iModelB, nBoot);
                    end
            end
        end
    end
end
