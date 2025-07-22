% filepath: /Users/sx712/Desktop/GitHub_local/PF_RC/Codes/shell_all_MC.m
% This script generates all combinations for OOD_MC_current and prints the command for each.
% You can adapt this to call a MATLAB function directly if desired.

flag_PatchMode = 2;
flag_cutMapping = 0;
flag_mirrorMapping = 1;
MCmode = 1;  % 1=10-fold CV, 2=LOOCV, 3=IC (AIC, AICc, BIC)

fprintf('Submitting jobs with: flag_PatchMode=%d, flag_cutMapping=%d, flag_mirrorMapping=%d, MCmode=%d\n', ...
    flag_PatchMode, flag_cutMapping, flag_mirrorMapping, MCmode);

for iSubj = 1:15
    for iiLoc_all_all = 1:3
        fprintf('Submitting: iSubj=%d, iiLoc_all_all=%d\n', iSubj, iiLoc_all_all);

        % ORI feature
        ifeature = 1; % ORI
        for ifamily = [1, 8] % 1=Gaussian, 8=DoG
            fprintf('  ifeature=%d (ORI), ifamily=%d\n', ifeature, ifamily);
            % Call your MATLAB function here, e.g.:
            OOD_MC_current(iSubj, iiLoc_all_all, flag_PatchMode, flag_cutMapping, flag_mirrorMapping, ifeature, ifamily, MCmode);
        end

        % SF feature
        ifeature = 2; % SF
        for ifamily = [2, 3] % 2=log parabola, 3=truncated log parabola
            fprintf('  ifeature=%d (SF), ifamily=%d\n', ifeature, ifamily);
            % Call your MATLAB function here, e.g.:
            OOD_MC_current(iSubj, iiLoc_all_all, flag_PatchMode, flag_cutMapping, flag_mirrorMapping, ifeature, ifamily, MCmode);
        end
    end
end