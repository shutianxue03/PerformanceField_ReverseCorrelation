
flag_stage = 1;
nBoot = 100;

for isubj = 1:15
    for iLocComb = [1,8,6,7,5,3]
        for iModelA = 1:2
            switch flag_stage
                case 1
                    OOD_NOM_Trialwise_compIV_temp(isubj, iLocComb, iModelA, nBoot);
                case 2
                    for iModelB = 1
                        OOD_NOM_Trialwise_fitNOM(isubj, iLocComb, iModelA, iModelB, nBoot);
                    end
            end
        end
    end
end
