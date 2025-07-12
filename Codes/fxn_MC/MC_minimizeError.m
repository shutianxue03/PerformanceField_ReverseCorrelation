

function err = MC_minimizeError(modelInd, params_fit_, xData, yData, model, fitMode, nparams_fullModel, ntrials_allSubj)
% fitMode: 1=SSE, 2=MLE

[nsubj, nLocTrue, nfilters] = size(yData);

% weight error terms by trial number
if nsubj>1, ntrialsProp = ntrials_allSubj'/sum(ntrials_allSubj); else, ntrialsProp = 1; end

ndata = nsubj * nLocTrue * nfilters;

% calculate the errors
err_ = nan(nsubj, nLocTrue);
for isubj = 1:nsubj
    
    yData_perSubj = squeeze(yData(isubj, :, :));
    
    % make predictions
    y_pred = RC_predKernel_MC(modelInd, xData, yData_perSubj, nparams_fullModel, params_fit_, model);
    
    % calculate error
    for iLoc = 1:nLocTrue
        switch fitMode
            case 1 % SSE
                err_(isubj, iLoc) = sum((yData_perSubj(iLoc, :) - y_pred(iLoc, :)).^2); 
            case 2 % MLE
                err_(isubj, iLoc) = -sum(log(normpdf(yData_perSubj(iLoc, :) , y_pred(iLoc, :), 1))); % MLE
        end
    end
    
end
err_ = ntrialsProp*err_;

if fitMode == 1
    err = sum(err_(:))/ndata;
else
    err = sum(err_(:));
end


