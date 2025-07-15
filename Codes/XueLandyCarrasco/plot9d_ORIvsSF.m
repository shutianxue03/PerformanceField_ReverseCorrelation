
switch iLocComb_all(1)
    case 6,
        nameAsym = 'HVA';
        load(sprintf('%s/%s/n%d_B1000_L67_%d_%d_m_ORI%dSF%d.mat', nameFolderCompile_RC, fileName, nsubj, nORI, nSF, ifamily_perF))
    case 5,
        nameAsym = 'VMA';
        load(sprintf('%s/%s/n%d_B1000_L53_%d_%d_m_ORI%dSF%d.mat', nameFolderCompile_RC, fileName, nsubj, nORI, nSF, ifamily_perF))
end

ifamily1 = ifamily_perF(1);
ifamily2 = ifamily_perF(1);
namesTunCs1 = namesTunC_unit_perF{ifamily1, paramMode};
namesTunCs2 = namesTunC_unit_perF{ifamily2, paramMode};
nTunCs_full1 = length(namesTunCs1);
nTunCs_full2 = length(namesTunCs1);
ifeature=2;
for iTunC = 2:nTunCs_full1
    
    x_allSubj = squeeze(margTuningC_ORI_allSubj(:, :, :, itype, iTunC));
    y_allSubj = squeeze(margTuningC_SF_allSubj(:, :, :, itype, iTunC));  % nsubj x nB x nLoc x ntypes x nTunC
    
    asymX_allSubj = (x_allSubj(:, :, 1)-x_allSubj(:, :, 2))./(x_allSubj(:, :, 1)+x_allSubj(:, :, 2));
    asymX_med_allSubj = getCI(asymX_allSubj, 1, 2);
    
    asymY_allSubj = (y_allSubj(:, :, 1)-y_allSubj(:, :, 2))./(y_allSubj(:, :, 1)+y_allSubj(:, :, 2));
    asymY_med_allSubj = getCI(asymY_allSubj, 1, 2);
    
    switch iLocComb_all(1)
        
        case 6 % HVA
            switch ifamily1
                case 1, ticks_allTunC_lb = -[100, 50, 50, 100]; ticks_allTunC_ub = [100, 50, 50, 100];
                case 8, ticks_allTunC_lb = -[60, 30, 45, 160, 20, 220]; ticks_allTunC_ub = [100, 40, 55, 160, 40, 220];
                case 2, ticks_allTunC_lb = -[50, 60, 50, 145]; ticks_allTunC_ub = [50, 60, 70, 265];
            end
        case 5 % VMA
            switch ifamily1
                case 1, ticks_allTunC_lb = -[100, 50, 50, 50]; ticks_allTunC_ub = [100, 50, 50, 50];
                case 8, ticks_allTunC_lb = -[100, 40, 300, 320, 50, 320]; ticks_allTunC_ub = [100, 44, 300, 320, 50, 320];
                case 2, ticks_allTunC_lb = -[50, 60, 60, 250]; ticks_allTunC_ub = [50, 60, 60, 250];
                    %                                 case 3, ticks_allTunC = [30, 50, 60, 110, 110];
            end
    end
    ticks = linspace(ticks_allTunC_lb(iTunC), ticks_allTunC_ub(iTunC), 5);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    text_title = sprintf('ORI vs. SF [%s]: %s', nameAsym, namesTunC_unit_perF{ifamily_perF(ifeature), 2}{iTunC});
    flag_sig = basicFxn_drawCorrAsym(asymX_med_allSubj*100, asymY_med_allSubj*100, ...
        ticks, ticks, ticks, ticks, text_title, markers_allSubj);
    saveas(gcf, sprintf('%s/L%d%d_TunC%d.jpg', folderNameFig_ORISF, iLocComb_all, iTunC))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end % iTunC
