function OOD_NOM_trialWise(isubj, iLocComb, iModelA, ni)

% generate IVs
fprintf('\n\n=====================\nGenerating IVs\n=====================\n\n')
nameFolder_NOM0 = 'Data_NOM_trialWise';
OOD_NOM_trialWise_beforeEst(isubj, iLocComb, iModelA, ni, nameFolder_NOM0);

% fit different models to IVs and make predictions
iModelB_all = 1:5; % 1=only constant noise; 2=only induced noise;3=no noise ()only lapse rate and criterion); 4=constant noise + criterion; 5=induced noise +criterion
for iModelB=iModelB_all
    fprintf('\n\n=====================\nModelB #%d\n=====================\n\n', iModelB)
    OOD_NOM_trialWise_Est(isubj, iLocComb, iModelA, iModelB, ni, nameFolder_NOM0)
end
end