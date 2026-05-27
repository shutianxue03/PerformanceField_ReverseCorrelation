

function deviance = MC_getDev_IDVD(modelInd, params_fit_, xData, yData, ifamily, fitMode)

[nlines, nfilters] = size(yData);
nparams_full = length(modelInd);

ndata = nlines * nfilters;

% empty container
deviance_ = nan(1, nlines);

% make predictions
y_pred = MC_predKernel(modelInd, xData, nlines, nparams_full, params_fit_, ifamily);

% calculate the deviance
for iLoc = 1:nlines
    switch fitMode
        case 1 % LS method
            deviance_(iLoc) = sum((yData(iLoc, :) - y_pred(iLoc, :)).^2); % no weight
%             if find(ifamily==[1,8])
%                 weights = exp(-(xData/1).^2);
%                 weights = weights/sum(weights);
%                 deviance_(iLoc) = sum((yData(iLoc, :) - y_pred(iLoc, :)) .^2.* weights); % weighted average
%             end
        case 2 % ML method
            deviance_(iLoc) = -sum(log(normpdf(yData(iLoc, :) , y_pred(iLoc, :), 1)));
    end
end

if fitMode == 1
    deviance = sum(deviance_(:))/ndata;
else
    deviance = sum(deviance_(:));
end
% deviance = sum(deviance_(:));


