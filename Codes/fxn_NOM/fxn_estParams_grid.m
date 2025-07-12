
% the OOD file with the same name
function nLL_3D = fxn_estParams_grid(sim, truth, fitMode, param_all_all, errorComp)

param1_all = param_all_all{1}; n1 = length(param1_all);
param2_all = param_all_all{2}; n2 = length(param2_all);
param3_all = param_all_all{3}; n3 = length(param3_all);

clc

fileName = sprintf('fxn_model/error_3D_%d_%d_%d.mat', n1, n2, n3);
fileDir = dir(fileName);

if isempty(fileDir)
    nLL_3D = nan(n1, n2, n3);
    for ip1 = 1:n1 % do not change to combvec, parfor does not accept that
        for ip2 = 1:n2
            parfor ip3 = 1:n3
                nLL_3D(ip1, ip2, ip3) = fxn_getError(fitMode, [param1_all(ip1), param2_all(ip2), param3_all(ip3)], sim, truth, errorComp);
            end
        end
        fprintf('ip1 = %d\n', ip1)
    end
    fprintf('DONE\n')
    
    save(fileName, 'nLL_3D')
    
else
    load(fileName)
    fprintf('LOADED\n')
end


