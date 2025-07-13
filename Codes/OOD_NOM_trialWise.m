function OOD_NOM_trialWise(isubj, iLocComb, iModelA, ni)

% OOD_NOM_trialWise: This function generates IVs and fits different models to them.
% Inputs:
%      isubj: index of subj
%      iLocComb: 1=Fovea, 8=periF(6 deg ecc), 6=HM, 7=VM, 5=LVM, 3=UVM
%      iModelA: 1=core model, 2=permuted template, 3=use IO template
%      ni: number of iterations

% generate IVs
fprintf('\n\n=====================\nGenerating IVs\n=====================\n\n')
nameFolder_Data = '/Volumes/purplab/EXPERIMENTS/1_Current_Experiments/Shutian_server/PF_RC/Data';
nameFolder_NOM0 = 'Data_NOM_trialWise';
OOD_NOM_trialWise_beforeEst(isubj, iLocComb, iModelA, ni, nameFolder_Data, nameFolder_NOM0);

% fit different models to IVs and make predictions
iModelB_all = 1:5; % 1=only constant noise; 2=only induced noise;3=no noise ()only lapse rate and criterion); 4=constant noise + criterion; 5=induced noise +criterion
for iModelB=iModelB_all
    fprintf('\n\n=====================\nModelB #%d\n=====================\n\n', iModelB)
    OOD_NOM_trialWise_Est(isubj, iLocComb, iModelA, iModelB, ni, nameFolder_Data, nameFolder_NOM0)
end
end