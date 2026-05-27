
function deviance = MC_getDev_CV(paramInd, params_est, x, y, ifamily, fitMode)
% MC_getDev_CV2 exists
% SX_normPDF = @(x,mu,sigma) exp(-(x-mu).^2/(2*sigma^2));

%%
[nlinesY, nfilters] = size(y);
% nparams_full = length(paramInd);
nparams_full = length(paramInd);
ndata = nlinesY * nfilters; % after reorganizing the matrix (y) to a vector

%% empty container
deviance = nan(1, nlinesY);

%% predict
y_pred = MC_predKernel(paramInd, x, nlinesY, nparams_full, params_est, ifamily);

% simple & quick plot
% figure('Position',[2e3 5e2 2e2 2e2]), hold on, 
% plot(x(1,:), y(1, :), 'or'), plot(x(1,:), y_pred(1, :), 'r-')
% plot(x(2,:), y(2, :), 'ob'), plot(x(2,:), y_pred(2, :), 'b-')

% calculate the deviance for each line
for iline = 1:nlinesY
    switch fitMode
        case 1 % LS method
            deviance(iline) = sum((y(iline, :) - y_pred(iline, :)).^2);
        case 2 % ML method
            deviance(iline) = -sum(log(normpdf(y(iline, :), y_pred(iline, :), 1)));
    end
end

%% get the combined deviance
if fitMode == 1
    deviance = sum(deviance(:))/ndata;
else
    deviance = sum(deviance(:));
end


