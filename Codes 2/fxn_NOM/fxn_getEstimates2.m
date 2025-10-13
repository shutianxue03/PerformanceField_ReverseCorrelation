
% assume that the IV distribution is PRS-ABS 

function pred = fxn_getEstimates2(params_est, IV_PRS, IV_ABS, ntrials_sim, normIV, plotFlag)
% get nLoc
nLoc = length(IV_PRS);
% normalize (all data centered at 0 and have SD=1)
for iLoc = 1:nLoc
    
    IV_diff = IV_PRS{iLoc} - IV_ABS{iLoc};
        IV_diff_mean = mean(IV_diff(:));
    IV_diff_std = std(IV_diff(:));
    IV_diff_norm = (IV_diff-IV_diff_mean)/IV_diff_std;
%     IV = [IV_PRS{iLoc}; IV_ABS{iLoc}];
%     IV_mean = mean(IV(:));
%     IV_std = std(IV(:));
%     IV_PRS_norm = (IV_PRS{iLoc} - IV_mean)/IV_std;
%     IV_ABS_norm = (IV_ABS{iLoc} - IV_mean)/IV_std;
   
    
    % calculate the mean/std of the IV distribution
    mean_IV_PRS_norm = mean(IV_PRS_norm);
    std_IV_PRS_norm = std(IV_PRS_norm);
    mean_IV_ABS_norm = mean(IV_ABS_norm);
    std_IV_ABS_norm = std(IV_ABS_norm);
    
    % plot histogram
    if plotFlag
        %     figure('Position', [0 0 1800 200])
        ModelPlot_hist(1, 1, IV_PRS_norm, IV_ABS_norm, params_est(iLoc, :))
    end
    
    % get metrics by simulation
    pred_perLoc = fxn_SimPred(ntrials_sim, mean_IV_PRS_norm, std_IV_PRS_norm, mean_IV_ABS_norm, std_IV_ABS_norm, params_est);
    
    x = linspace(0, .5, 1e4);
    thresh = params_est(3);
    pdf_PRS_ = normpdf(x-thresh, mean_IV_PRS_norm, std_IV_PRS_norm); pdf_PRS = pdf_PRS_/sum(pdf_PRS_);
    pdf_ABS_ = normpdf(x-thresh, mean_IV_ABS_norm, std_IV_ABS_norm);pdf_ABS = pdf_ABS_/sum(pdf_ABS_);
    cdf_PRS = normcdf(x-thresh, mean_IV_PRS_norm, std_IV_PRS_norm);
    cdf_ABS = normcdf(x-thresh, mean_IV_ABS_norm, std_IV_ABS_norm);
    pA_PRS_math = sum((cdf_PRS.^2 + (1-cdf_PRS).^2) .* pdf_PRS);
    pA_ABS_math =  sum(((1-cdf_ABS).^2 + cdf_ABS.^2) .* pdf_ABS);
    pA_math = mean([pA_PRS_math, pA_ABS_math]);
    
    % assign to each loc
    pred.IV_PRS{iLoc} = pred_perLoc.IV_PRS;
    pred.IV_ABS{iLoc} = pred_perLoc.IV_PRS;
    pred.dprime(iLoc) = pred_perLoc.dprime;
    pred.criterion(iLoc) = pred_perLoc.criterion;
    pred.pC(iLoc) = pred_perLoc.pC;
    pred.pHit(iLoc) = pred_perLoc.pHit;
    pred.pFA(iLoc) = pred_perLoc.pFA;
%     pred.pA(iLoc) = pred_perLoc.pA;
%     pred.pA_PRS(iLoc) = pred_perLoc.pA_PRS;
%     pred.pA_ABS(iLoc) = pred_perLoc.pA_ABS;
pred.pA(iLoc)  = pA_math;
pred.pA_PRS(iLoc) = pA_PRS_math;
pred.pA_ABS(iLoc) = pA_ABS_math;
end

end