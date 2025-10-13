
function [IV_PRS, IV_ABS, maxPRS_allT, maxABS_allT] = fxn_getIV_v2(iModelA, e2D, iPRS, convolveType, IVType, template_true, ORI_bound)


%% extract energy
e_perT_PRS = squeeze(e2D(boolean(iPRS), :, :));
e_perT_ABS = squeeze(e2D(boolean(1-iPRS), :, :));

%% get IV and max
[maxPRS_allT, IV_PRS] = fxn_IV_perPRSABS(e_perT_PRS, iModelA, convolveType, IVType, template_true, ORI_bound);
[maxABS_allT, IV_ABS] = fxn_IV_perPRSABS(e_perT_ABS, iModelA, convolveType, IVType, template_true, ORI_bound);

%%
function [max_allT, IV] = fxn_IV_perPRSABS(e_perT, iModelA, convolveType, IVType, template_true, ORI_bound)
[ntrials, nORI, nSF] = size(e_perT);
IV = nan(ntrials, 1);
max_allT = nan(ntrials, 2);

% randomize pixels of template
% (though there is only ONE template for all trials, still need to randomize it differently on each trial)
parfor itrial = 1:ntrials
    try
        if (iModelA == 2) || (iModelA == 4) % randomize template on each trial
            rng('shuffle')
            indRand = randperm(nORI*nSF);
            template_true_v = template_true(:);
            template_rand = template_true_v(indRand);
            template = reshape(template_rand, nORI, nSF);
        else
            template = template_true;
        end
        
    catch, fprintf('Alert: problem with TEMPLATE (line 23-31)\n')
    end
    
    try
        e = squeeze(e_perT(itrial, :, :));
    catch, fprintf('Alert: problem with ENERGY (line 37)\n')
    end
    
    try
        switch convolveType
            case 1, tempCovE = template.*e; % dot multiply
            case 2 , tempCovE = conv2(template, e, 'same'); % convolve
        end
    catch, fprintf('Alert: problem with integration  (line 41-44)\n')
    end
    
    % get the channel at which max value is found
    try
        plotFlag = 0;
%         [max_, maxORI, maxSF_log] = getIndMax(tempCovE, ORI_bound, plotFlag);
        maxORI=nan;maxSF_log=nan;
        max_allT(itrial, :) = [maxORI, maxSF_log];
    catch, fprintf('Alert: problem with getting max (line 50-52)\n')
    end
    
    %     quickPlot_tempCovE, waitforbuttonpress
    
    try
        % decide the format of IV
        switch IVType
            case 1, IV(itrial) = sum(tempCovE(:));
            case 2, IV(itrial) = max_;
            case 3, IV(itrial) = max_/(sum(tempCovE(:))); if sum(tempCovE(:)) == 0, error('ALERT: The suppressive drive is 0.'); end% SS stands for the estimated semi-saturation factor
        end
    catch, fprintf('Alert: problem with IV (line 59-63)\n')
    end
end % end of itrial