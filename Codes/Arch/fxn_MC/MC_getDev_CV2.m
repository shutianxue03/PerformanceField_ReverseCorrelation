

function deviance = MC_getDev_CV2(params_est, x, y, ifamily, fitMode)

% fit each loc separately
% SX_normPDF = @(x,mu,sigma) exp(-(x-mu).^2/(2*sigma^2));

%%
nfilters = length(y);

%% predict
y_pred = predSFkernel(x, ifamily, params_est, 0);
%% simple & quick plot
% figure('Position',[2e3 5e2 2e2 2e2]), hold on,
% plot(x, y, 'or'), plot(x(1,:), y_pred(1, :), 'r-')

%% calculate the deviance for each line
switch fitMode
    case 1 % LS method
        deviance = sum((y(:) - y_pred(:)).^2);
    case 2 % ML method
        deviance = -sum(log(normpdf(y(:), y_pred(:), 1)));
end

%% get the combined deviance
if fitMode == 1
    deviance = sum(deviance(:))/nfilters;
else
    deviance = sum(deviance(:));
end

end


