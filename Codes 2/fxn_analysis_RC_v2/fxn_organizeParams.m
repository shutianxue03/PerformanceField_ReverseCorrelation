
function params_fit = fxn_organizeParams(printFlag, namesParams, nParams_full, nLocTrue, modelInd, params_fit_)
params_fit = cell(1,nParams_full);
paramsLen = nLocTrue .^ modelInd; %

if printFlag, fprintf('%s\n', num2str(modelInd)), end

for iParam = 1:nParams_full
    % decide the start
    if iParam==1, iStart = 1;
    else, iStart = sum(paramsLen(1:iParam-1))+1;
    end
    
    % decide the end
    if iParam == nParams_full, iEnd =  sum(paramsLen);
    else, iEnd = sum(paramsLen) - sum(paramsLen(iParam+1:end));
    end
    
    params_fit{iParam} = params_fit_(iStart:iEnd);
    
    if printFlag
        fprintf('%s: %s\n', namesParams{iParam}, num2str(params_fit{iParam}))
    end
end

if printFlag, fprintf('\n'), end

