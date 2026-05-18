noiseCST = .2;
gaborCST = .3;
nTrials = 1e4;
nIter = 50;
cSDT_true = 0;
lambda_whiten = 0;
% nBasisORI = 6;
% basisWidthORI = 0.9;
% nBasisSF = 6;
% basisWidthSF = 0.6;

for iModelB_sim = 6
    for nBasisORI = [3,9]
         for basisWidthORI = [.5, 1.1]
             for nBasisSF =  [3,9]
                 for basisWidthSF = [.5, 1.1]
                     disp(['iModelB_sim: ' num2str(iModelB_sim) ', nBasisORI: ' num2str(nBasisORI) ', basisWidthORI: ' num2str(basisWidthORI) ', nBasisSF: ' num2str(nBasisSF) ', basisWidthSF: ' num2str(basisWidthSF)])
                     OOD_sim(noiseCST, gaborCST, nTrials, Nmul_true, Nadd_true, Nshared_true, cSDT_true, lambda_whiten, iModelB_sim, nIter, nBasisORI, basisWidthORI, nBasisSF, basisWidthSF)
                 end
             end
         end
    end
end
