
function y_pred = MC_predKernel(paramInd, x, nlinesY, nparams_full, params_est, ifamily)

% paramInd: a vector indicating whether the param is free to vary (1) or shared across locations (0)
% nlinesY: number of lines in y
% nparams_full: number of parameters considered
% params_est: contains both fixed and shared parameter of the two (or multiple) locations

nfilters = size(x, 2);

%% reorganize params
printFlag = 0;
params_est_re = fxn_organizeParams(printFlag, {}, nparams_full, 2, paramInd, params_est);

%%
y_pred = nan(nlinesY, nfilters);

for iline = 1:nlinesY
    params = nan(1,nparams_full);
    
    for iparam = 1:nparams_full
        % assign the chosen parameters to each loc
        if paramInd(iparam) == 1
            params(iparam) = params_est_re{iparam}(iline);
        else
            params(iparam) = params_est_re{iparam};
        end
    end
    
%     if iparam == 5, params(5) = 1-params(5); end
    if  size(x,1) == 1
        y_pred(iline, :) = predSFkernel(x, ifamily, params, 0);
    else
        y_pred(iline, :) = predSFkernel(x(iline, :), ifamily, params, 0);
    end
end


