clc, close all

subjList =            {'YK', 'SP', 'SX', 'LS', 'RE', 'MD',         'HL', 'FH', 'HA',  'DT', 'CS', 'DU',        'SR'}; % no AS or RC, n=13
nblocks_allSubj = [200, 240, 210, 240, 220, 220,        220, 205, 210, 205, 205, 205,      195];
nsubj = length(subjList);

nLoc = 4; % left HM, UVM, right HM, LVM
pC_ave_allSubj = nan(nsubj, nLoc);
cst_log_ave_allSubj = pC_ave_allSubj;
CS_ave_allSubj = pC_ave_allSubj;
pC_std_allSubj = pC_ave_allSubj;
cst_log_std_allSubj = pC_ave_allSubj;
CS_std_allSubj = pC_ave_allSubj;

for isubj = 1:nsubj
    subjName = subjList{isubj};
    nblocks = nblocks_allSubj(isubj);
    % load nehavMeas
    
    load(sprintf('Data_OOD/%s%d/%s_behavMeas.mat', subjName, nblocks, subjName), 'pC3_perSess_perLoc', 'cst_perSess_perLoc');
    % pC3_perSess_perLoc: nSess x nLoc x 3
    % cst_perSess_perLoc: nSess x nLoc
    pC_allSess = squeeze(pC3_perSess_perLoc(:, 2:5, 1));
    cst_ln_allSess = squeeze(cst_perSess_perLoc(:, 2:5, 1));
    
    cst_log_allSess = log10(cst_ln_allSess);
    CS_allSess = 1./cst_ln_allSess;
    
    pC_ave_allSubj(isubj, :) = mean(pC_allSess);
    cst_log_ave_allSubj(isubj, :) = mean(cst_log_allSess);
    CS_ave_allSubj(isubj, :) = mean(CS_allSess);
    pC_std_allSubj(isubj, :) = std(pC_allSess);
    cst_log_std_allSubj(isubj, :) = std(cst_log_allSess);
    CS_std_allSubj(isubj, :) = std(CS_allSess);
    
    figure
    subplot(3,1,1), plot(pC_allSess), yline(.7, 'k-'); ylim([.5, 1])
    subplot(3,1,2), plot(cst_log_allSess)
    subplot(3,1,3), plot(CS_allSess)
    sgtitle(subjName)
end
%%
clc
round(mean(pC_std_allSubj)*100, 1)
round(std(pC_std_allSubj)/sqrt(nsubj)*100, 1)

round(mean(CS_std_allSubj)*100, 1)
round(std(CS_std_allSubj)/sqrt(nsubj)*100, 1)

indLoc = repmat(1:nLoc, nsubj, 1) ; print_nANOVA({'Loc'}, pC_std_allSubj(:), indLoc(:), nsubj)
indLoc = repmat(1:nLoc, nsubj, 1) ; print_nANOVA({'Loc'}, CS_std_allSubj(:), indLoc(:), nsubj)