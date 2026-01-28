
% flag_stage = 1;
nIter = 10;
nJob = 2;

for flag_stage=1:2
    for iJob=1:nJob
        for isubj = 1:2
            for iLocComb = 1:5
                for iModelA = 1:2
                    switch flag_stage
                        case 1
                            OOD_NOM_Trialwise_compIV(isubj, iLocComb, iModelA, nIter, nJob, iJob);
                        case 2
                            for iModelB = 1:2
                                OOD_NOM_Trialwise_fitNOM(isubj, iLocComb, iModelA, iModelB, nIter, nJob, iJob);
                            end
                    end
                end
            end
        end
    end
end