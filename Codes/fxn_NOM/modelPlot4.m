
sz_fig = [350 500];

for iiLoc = 1:nLoc
    iLocComb = iLocComb_all(iiLoc);
    
    switch iLocComb
        case 6, yticks_ = linspace(0, 120, 5);
        case 7, yticks_ = linspace(0, 240, 5);
        case 5, yticks_ = linspace(0, 240, 5);
        case 3, yticks_ = linspace(0, 240, 5);
            %
        case 6, yticks_ = linspace(0, 4, 5);
        case 5, yticks_ = linspace(0, 12, 5);
        case 3, yticks_ = linspace(0, 6, 5);
    end
    
    %     yticks_ = linspace(0, 1.2e3, 5);
    for iIC = 3%1:nIC % all IC shows the same value
        
        % get delta IC
        switch flag_compareAcrossA_OR_B
            case 1
                %                 [IC_med, ~, ~, IC_CI_neg, IC_CI_pos] = getCI(NOM_nLL_allSubj(indSubj, iiLoc, :, iModelB, :), 1, 5); % compare across A types: core, random energy, rand. temp, rand. both
                [IC_med, ~, ~, IC_CI_neg, IC_CI_pos] = getCI(NOM_IC_allSubj(:, iiLoc, :, iModelB, :), 1, 5); % compare across A types: core, random energy, rand. temp, rand. both
            case 2
                iModelA = 1; [IC_med, ~, ~, IC_CI_neg, IC_CI_pos] = getCI(NOM_nLL_allSubj(indSubj, iiLoc, iModelA, :, :), 1, 5); % compare across B types: full, constantN only, inducedNonly
        end
        
        if size(IC_med, 2)==1, IC_med_delta = IC_med;
        else,
            IC_med_delta = IC_med - min(IC_med, [], 2);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        switch flag_compareAcrossA_OR_B
            case 1
                basicFxn_drawBars(IC_med_delta, [], repmat(colors_comb(iLocComb, :), nModelsA, 1), namesModelA(1:nModelsA), ...
                    yticks_, yticks_, flag_plotIDVD, 0, sprintf('L%d-%s', iLocComb, namesIC{iIC}), 1, sz_fig);
                %                     yticks_, yticks_, flag_plotIDVD, 0, sprintf('L%d-%s', iLocComb, namesIC{iIC}), 1, sz_fig);
                
            case 2
                basicFxn_drawBars(IC_med_delta, [], repmat(colors_comb(iLocComb, :), nModelsB, 1), namesModelB(1:nModelsB), ...
                    yticks_, yticks_, flag_plotIDVD, 0, sprintf('L%d-%s', iLocComb, namesIC{iIC}), 1, sz_fig);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        subjList(IC_med_delta(:, 1)~=0)
        % save
        folderName = sprintf('%s/temp%d/modelVariations/', nameFolderFig, templateType);
        folderDir = dir(folderName); if isempty(folderDir), mkdir(folderName), end
        saveas(gcf, sprintf('%sn%d_L%d_%s.jpg', folderName, nsubj, iLocComb, namesIC{iIC}))
        
    end % iIC
    
end % iiLoc

