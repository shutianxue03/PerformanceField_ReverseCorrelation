function OOD_NOM_Trialwise(isubj, iLocComb, iModelA,nIterations)

% OOD_NOM_trialWise: This function generates IVs and fits different models to them.
% Inputs:
%      isubj: index of subj
%      iLocComb: 1=Fovea, 8=periF(6 deg ecc), 6=HM, 7=VM, 5=LVM, 3=UVM
%      iModelA: 1=core model, 2=permuted template, 3=use IO template
%      ni: number of iterations

%% Generate IVs
fprintf('\n\n=====================\nGenerating IVs\n=====================\n\n')

% The settings below can be placed in SX_RC1_setting
IVType = 1;            % 1=sum of the dot product/convolution; 2=max; 3=normalized
templateType = 1; % (1) raw (2) reconstructed kernel (3) mirrored template
flag_fminconORbads = 2; % 1=use fmincon when fitting NOM to data, faster; 2=bads, slower but better

%------------------------------%
OOD_NOM_Trialwise_beforeEst(isubj, iLocComb, iModelA, IVType, templateType, nIterations);
%------------------------------%

%% Fit different models to IVs and make predictions
iModelB_all = 5; % 4=only SDadd; 5=only Nmul
for iModelB=iModelB_all
    fprintf('\n\n=====================\nModelB #%d\n=====================\n\n', iModelB)
    %------------------------------%
    OOD_NOM_Trialwise_Est(isubj, iLocComb, iModelA, iModelB, nIterations, flag_fminconORbads)
    %------------------------------%
end

