
colors = colors_comb(iLocComb_all, :);

% for iModelB = 1:nModelsB
clear *med_allSubj *neg_allSubj *pos_allSubj
namesParamsNOM = namesParamsModel_all{iModelB};
nNoise = length(namesParamsNOM);

for iModelA = 1%indModelsA%:nModelsA % core model only
    
    % empty containers
    NOM_data_med_allSubj = nan(nsubj, nLoc, 8);
    NOM_data_neg_allSubj = NOM_data_med_allSubj;
    NOM_data_pos_allSubj = NOM_data_med_allSubj;
    NOM_pred_med_allSubj = NOM_data_med_allSubj;
    NOM_pred_neg_allSubj = NOM_data_med_allSubj;
    NOM_pred_pos_allSubj = NOM_data_med_allSubj;
    NOM_params_med_allSubj = nan(nsubj, nLoc, nNoise);
    NOM_params_neg_allSubj = NOM_params_med_allSubj;
    NOM_params_pos_allSubj = NOM_params_med_allSubj;
    
    for isubj = 1:nsubj
        for iiLoc = 1:nLoc
            [data_med, ~, ~, data_neg, data_pos] = getCI(NOM_data_metrics_allSubj(isubj, iiLoc, iModelA, iModelB, :, :), 1, 5);
            [pred_med, ~, ~, pred_neg, pred_pos] = getCI(NOM_pred_metrics_allSubj(isubj, iiLoc, iModelA, iModelB, :, :), 1, 5);
            [params_med, ~, ~, params_neg, params_pos] = getCI(NOM_params_est_allSubj(isubj, iiLoc, iModelA, iModelB, :, 1:nNoise), 1, 5);
            
            assert(length(params_med) == nNoise)
            NOM_data_med_allSubj(isubj, iiLoc, :) = data_med; NOM_data_neg_allSubj(isubj, iiLoc, :) = data_neg; NOM_data_pos_allSubj(isubj, iiLoc, :) = data_pos;
            NOM_pred_med_allSubj(isubj, iiLoc, :) = pred_med; NOM_pred_neg_allSubj(isubj, iiLoc, :) = pred_neg; NOM_pred_pos_allSubj(isubj, iiLoc, :) = pred_pos;
            NOM_params_med_allSubj(isubj, iiLoc, :) = params_med; NOM_params_neg_allSubj(isubj, iiLoc, :) = params_neg; NOM_params_pos_allSubj(isubj, iiLoc, :) = params_pos;
        end % iiLoc
    end % isubj
    
    %% Fig 1: metrics data vs. pred
    close all
    %----------%
    modelPlot1
    %----------%
    
    %% Fig 2: pA vs. params
    close all
    %----------%
    modelPlot2 % main figure presents the zero-meaned values and inset is the raw data
    %----------%
    
    %% Fig 3. compare params across loc
    close all
    flag_plotIDVD= 0;
    flag_pairwiseComp = 1;
    %----------%
    modelPlot3
    %----------%
    
    %     end % iModelB
end % iModelA

