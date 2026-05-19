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

Nmul_true=.5;
Nadd_true=5;
Nshared_true=5;

for iModelB_sim = 6
    for nBasisORI = [5,7]
         for basisWidthORI = [.7, .9]
             for nBasisSF =  [5,7]
                 for basisWidthSF = [.7, .9]
                     disp(['iModelB_sim: ' num2str(iModelB_sim) ', nBasisORI: ' num2str(nBasisORI) ', basisWidthORI: ' num2str(basisWidthORI) ', nBasisSF: ' num2str(nBasisSF) ', basisWidthSF: ' num2str(basisWidthSF)])
                     OOD_sim(noiseCST, gaborCST, nTrials, Nmul_true, Nadd_true, Nshared_true, cSDT_true, lambda_whiten, iModelB_sim, nIter, nBasisORI, basisWidthORI, nBasisSF, basisWidthSF)
                 end
             end
         end
    end
end
