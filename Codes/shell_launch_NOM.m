
flag_stage=1;
nIter=100;
nJob=10;

isubjList=1:15
iLocCombList=1:5
iModelAList=[1,2]
iModelBList=[1,2,3,4]
flag_whitenDV = 1;

for iJob=1:nJob
    for isubj = isubjList
        for iLocComb = iLocCombList
            for iModelA = iModelAList
                switch flag_stage
                    case 1
                        OOD_NOM_Trialwise_compIV(isubj, iLocComb, iModelA, nIter, nJob, iJob);
                    case 2
                        for iModelB = iModelBList
                            OOD_NOM_Trialwise_fitNOM(isubj, iLocComb, flag_whitenDV, iModelA, iModelB, nIter, nJob, iJob);
                        end
                end
            end
        end
    end
end
