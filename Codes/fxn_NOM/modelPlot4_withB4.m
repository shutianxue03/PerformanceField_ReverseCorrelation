
sz_fig = [250 430];

% get average of random models
% iModelA: 1=core, 2=rand T, 3=rand E, 4=rand both
% iModelB: 1=Core noisy model (with internal noise components); 4=Core noiseless model (no internal noise components)

% NOM_nLL_allSubj: nsubj x nLoc x nModels x 1 x ni
% NOM_IC_allSubj: nsubj x nLoc x nModels x 1 x ni x nIC
if iIC==0, NOM_GoF_allSubj = squeeze(NOM_nLL_allSubj);
else, NOM_GoF_allSubj = squeeze(NOM_IC_allSubj(:, :, : , :, :, iIC));
end
% NOM_GoF_allSubj: nsubj x nLoc x nModelA x ni
NOM_GoF_allSubj_core = squeeze(NOM_GoF_allSubj(:, :, 1 , :)); % nsubj x nLoc x ni
NOM_GoF_allSubj_rand_all = squeeze(NOM_GoF_allSubj(:, :, 2:4 , :)); % nsubj x nLoc x 3 x ni

% take the average of all three random models
NOM_GoF_allSubj_rand_selected = getCI(NOM_GoF_allSubj_rand_all, 2, 3); % nsubj x nLoc x ni
% or just one random model (1=template, 2=energy profile, 3=both)
% NOM_GoF_allSubj_rand_selected = squeeze(NOM_GoF_allSubj_rand_all(:, :, 3, :)); % nsubj x nLoc x ni

NOM_GoF_allSubj_B1_med = getCI(cat(4, NOM_GoF_allSubj_core, NOM_GoF_allSubj_rand_selected), 1, 3); % nsubj x nLoc x 2

% get the inefficient model
% nameFile_B4 = input(sprintf('The B1 file name is %s\nEnter the B4 file name: ', nameFileNOM));
nameFile_B4 = 'Data_compile/NOM/ORI29SF29/n13_n1000_Lall4_N_temp3_B4_nA1_conv1_IV1.mat';
load(nameFile_B4)
iModelA_B4 = 1;
iModelB_B4 = 4;
if iIC==0, NOM_GoF_allSubj_B4 = squeeze(NOM_nLL_allSubj(:, :, iModelA_B4, iModelB_B4, :)); % nsubj x nLoc x ni
else, NOM_GoF_allSubj_B4 = squeeze(NOM_IC_allSubj(:, :, iModelA_B4, iModelB_B4, :, iIC)); % nsubj x nLoc x nModels x 1 x ni x nIC
end

% NOM_GoF_allSubj_B4 = squeeze(NOM_GoF_allSubj(:, :, iModelA_B4, iModelB_B4, :, iIC)); % nsubj x nLoc x ni
NOM_GoF_allSubj_B4_med = getCI(NOM_GoF_allSubj_B4, 1, 3); % nsubj x nLoc

% core vs. noiseless model vs. rand model
NOM_GoF_allSubj_allM_med = cat(3, NOM_GoF_allSubj_B1_med(:, :, 1), NOM_GoF_allSubj_B4_med, NOM_GoF_allSubj_B1_med(:, :, 2));

% core vs. rand model
% NOM_GoF_allSubj_allM_med = NOM_GoF_allSubj_B1_med;

% core vs. noiseless model
% NOM_GoF_allSubj_allM_med = cat(3, NOM_GoF_allSubj_B1_med(:, :, 1), NOM_GoF_allSubj_B4_med);

%%
[nsubj, nLoc, nModelsAB] = size(NOM_GoF_allSubj_allM_med);

indLoc = nan(size(NOM_GoF_allSubj_allM_med));
indMod = indLoc;

for iLoc=1:nLoc
    indLoc(:, iLoc,:) = ones(nsubj, nModelsAB) * iLoc;
end

for iMod=1:nModelsAB
    indMod(:, :, iMod)=ones(nsubj, nLoc) * iMod;
end

text = print_nANOVA({'Loc', 'Model'}, NOM_GoF_allSubj_allM_med(:), {indLoc(:), indMod(:)}, nsubj)

%% PLOT
for iiLoc = 1:nLoc
    iLocComb = iLocComb_all(iiLoc);
    
    if iIC==0
        switch iLocComb
            case 6, yticks_ = linspace(0, 240, 5)/2;
            case 7, yticks_ = linspace(0, 240, 5)/2;
            case 5, yticks_ = linspace(0, 240, 5)/2;
            case 3, yticks_ = linspace(0, 240, 5)/2;
        end
    else
        switch iLocComb
            case 6, yticks_ = linspace(0, 120, 5);
            case 7, yticks_ = linspace(0, 240, 5);
            case 5, yticks_ = linspace(0, 240, 5);
            case 3, yticks_ = linspace(0, 240, 5);
        end
    end
    if flag_plotIDVD, yticks_ = yticks_*3; end
    GoF_med = squeeze(NOM_GoF_allSubj_allM_med(:, iiLoc, :));
    xticklabels_ =  {'CoreNoisy', 'CoreNoiseless', 'Random'}; assert(length(xticklabels_ ) == size(GoF_med, 2));
%     xticklabels_ =  {'Core',  'Random'}; assert(length(xticklabels_ ) == size(GoF_med, 2));
%     xticklabels_ =  {'Core',  'Inefficient'}; assert(length(xticklabels_ ) == size(GoF_med, 2));
    
    GoF_med = squeeze(NOM_GoF_allSubj_allM_med(:, iiLoc, :));
    GoF_med_delta = GoF_med - min(GoF_med, [], 2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if iIC==0
        basicFxn_drawBars(GoF_med_delta, [], repmat(colors_comb(iLocComb, :), nModelsA+1, 1), xticklabels_, ...
            yticks_, yticks_, flag_plotIDVD, 0, sprintf('L%d-%s', iLocComb, 'nLL'), 1, sz_fig, nPairs, str_tail);
        
    else
        basicFxn_drawBars(GoF_med_delta, [], repmat(colors_comb(iLocComb, :), nModelsA+1, 1), xticklabels_, ...
            yticks_, yticks_, flag_plotIDVD, 0, sprintf('L%d-%s', iLocComb, namesIC{iIC}), 1, sz_fig, nPairs, str_tail);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    subjList(GoF_med_delta(:, 1)>GoF_med_delta(:, 3))
    ind=GoF_med_delta(:, 1)>GoF_med_delta(:, 3); GoF_med_delta(ind, 1)-GoF_med_delta(ind, 3)
    % save
    folderName = sprintf('%s/temp%d/modelVariations/', nameFolderFig, templateType);
    folderDir = dir(folderName); if isempty(folderDir), mkdir(folderName), end
    if iIC==0
        saveas(gcf, sprintf('%sn%d_L%d_nLL.jpg', folderName, nsubj, iLocComb))
    else
        saveas(gcf, sprintf('%sn%d_L%d_%s.jpg', folderName, nsubj, iLocComb, namesIC{iIC}))
    end
    
end % iiLoc

