function fxn_plotPred(iLoc_all, nsubj, data, pred, params_est, suptitle_)
%  fxn_plotPred(nLoc, nsubj, data, pred, params_est, suptitle_)
% data structure
%       data: each field is nLoc x 1 (IDVD) or nsubj x nLoc (Group)
%       pred: each field is ni/nsubj x nLoc
%       params_est: ni/nsubj x nLoc x nparams

nparams_model = size(params_est, 3);

figure('Position', [0 0 1000 1000])
plotPred(1, nsubj, iLoc_all, data.dprime, pred.dprime, 'dprime'), legend({'Data' 'Prediction'}, 'Location', 'best')
plotPred(2, nsubj, iLoc_all, data.criterion, pred.criterion, 'criterion', [-1,1])
plotPred(4, nsubj, iLoc_all, data.pC, pred.pC, 'pC', [.5,1])
plotPred(5, nsubj, iLoc_all, data.pHit, pred.pHit, 'pHit', [.5,1])
plotPred(6, nsubj, iLoc_all, data.pFA, pred.pFA, 'pFA', [0, 1])
plotPred(7, nsubj, iLoc_all, data.pA, pred.pA, 'pA', [.5,1])
plotPred(8, nsubj, iLoc_all, data.pA_PRS, pred.pA_PRS, 'pA [PRS]', [.5,1])
plotPred(9, nsubj, iLoc_all, data.pA_ABS, pred.pA_ABS, 'pA [ABS]', [.5,1])
plotPred(10, nsubj, iLoc_all, zeros(size(data.dprime)), squeeze(params_est(:, :,1)), 'alpha [PRS]', [0, 1.5])
plotPred(11, nsubj, iLoc_all, zeros(size(data.dprime)), squeeze(params_est(:, :,2)), 'alpha [ABS]', [0, 1.5])
plotPred(12, nsubj, iLoc_all, zeros(size(data.dprime)), squeeze(params_est(:, :,3)), 'Threshold', [])

sgtitle(suptitle_)

