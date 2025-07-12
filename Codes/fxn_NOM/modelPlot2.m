
% Fig 2: pA vs. NOM params

flag_zeroMean = 1;
typeCorr = 'Pearson';
type_tail = 'both';

x_ticks = linspace(-.14, .14, 5);
y_ticks_all = [1.2, .16, nan]; % symmetrical

% get x
im_pA = 6;
x_med = NOM_data_med_allSubj(:, :, im_pA); % x: pA
% x_neg = NOM_data_neg_allSubj(:, :, im_pA); % x: pA
% x_pos= NOM_data_pos_allSubj(:, :, im_pA); % x: pA

for iNOMp = 1:(nNoise-1)%1:3

    y_med = NOM_params_med_allSubj(:, :, iNOMp);
%     y_neg = NOM_params_neg_allSubj(:, :, ip);
%     y_pos = NOM_params_pos_allSubj(:, :, ip);
    y_ticks = linspace(-y_ticks_all(iNOMp), y_ticks_all(iNOMp), 5);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    basicFxn_drawCorr(x_med, y_med, colors_comb(iLocComb_all, :), ...
        x_ticks, y_ticks, x_ticks, y_ticks, flag_zeroMean, typeCorr, type_tail, namesParamsNOM{iNOMp}, markers_allSubj);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % save
    folderName = sprintf('%s/temp%d/corr_pA_params/%s/', nameFolderFig, templateType, nameFileLoc);
    folderDir = dir(folderName); if isempty(folderDir), mkdir(folderName), end
    saveas(gcf, sprintf('%sn%d_NOMp%d_A%dB%d.jpg', folderName, nsubj, iNOMp, iModelA, iModelB))
    
end % ip


