
function [IV_PRS, IV_ABS, imax_allT] = fxn_getIV(ntrials, e2D, convolveType, IVType, template, ORI_bound)
%% empty containers
IV_PRS = nan(ntrials, 1);
IV_ABS = IV_PRS;
imax_allT = nan(ntrials, 2, 2); % ntrials x ORI/SF x PRS/ABS

%% calculate internal variable
for itrial = 1:ntrials
    e_perT_PRS = squeeze(e2D(itrial, :, :));
    e_perT_ABS = squeeze(e2D(itrial+ntrials, :, :));
    
    switch convolveType.i
        case 1 % dot multiply
            tempCovE_PRS = template.*e_perT_PRS;
            tempCovE_ABS = template.*e_perT_ABS;
        case 2 % convolve
            tempCovE_PRS = conv2(template, e_perT_PRS, 'same');
            tempCovE_ABS = conv2(template, e_perT_ABS, 'same');
    end
    
    %% get the channel at which max value is found    
%     if itrial == 1
%         plotFlag = 1;
%         fprintf('Max value is restricted with ORI = +/- 40 deg (manully set index!!)\n')
%     else, 
        plotFlag = 0;
%     end

    
    if plotFlag
        figure('Position', [2000 200 800 300])
        subplot(1,2,1), title('PRS')
    end
    [max_PRS, imax_ORI_PRS, imax_SF_PRS] = getIndMax(tempCovE_PRS, ORI_bound, plotFlag);
    
    if plotFlag, subplot(1,2,2), title('ABS'), end
    [max_ABS, imax_ORI_ABS, imax_SF_ABS] = getIndMax(tempCovE_ABS, ORI_bound, plotFlag);
    
    imax_allT(itrial, :, :) = [imax_ORI_PRS, imax_ORI_ABS; imax_SF_PRS, imax_SF_ABS];
    
    %% decide the format of IV
    switch IVType.i
        case 1
            IV_PRS(itrial) = sum(tempCovE_PRS(:));
            IV_ABS(itrial) = sum(tempCovE_ABS(:));
        case 2
            IV_PRS(itrial) = max_PRS;
            IV_ABS(itrial) = max_ABS;
        case 3 % normalization (Fernandez 2022)
            IV_PRS(itrial) = max_PRS/(sum(tempCovE_PRS(:)) + SS_PRS); % SS stands for semi-saturation
            IV_ABS(itrial) = max_ABS/(sum(tempCovE_ABS(:)) + SS_ABS);
    end
end % end of itrial
