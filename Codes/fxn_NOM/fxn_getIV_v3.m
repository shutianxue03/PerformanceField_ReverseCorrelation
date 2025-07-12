
function [IV, e3D_noisy, max_allT] = fxn_getIV_v3(iModelA, e3D, convolveType, IVType, template_true, ORI_bound)

[nTrials, nORI, nSF] = size(e3D);
IV = nan(nTrials, 1);
max_allT = nan(nTrials, 2);

% Add Gaussian noise to the template, controlled by sigma_template.
% This is done by multiplying a Gaussian matrix (sigma = 1) with 
% the energy profiles, then adding it to the integration between 
% the true template and the energy profile.
% e3D_noisy = e3D .* randn(size(e3D));

% CANNOT permutate pixels of template 
% to have the same IV for paired trials, the template cannot be randomized
% independently across trials (at least not for the pair of trials)

if iModelA == 2 % Permute pixels in the template
    rng('shuffle')
    indRand = randperm(nORI*nSF);
    template_true_v = template_true(:);
    template_rand = template_true_v(indRand);
    template = reshape(template_rand, nORI, nSF);
else
    template = template_true;
end

parfor itrial = 1:nTrials
    e = squeeze(e3D(itrial, :, :));
    
    switch convolveType
        case 1, tempCovE = template.*e; % cross-correlation
        case 2 , tempCovE = conv2(template, e, 'same'); % convolution
    end
    
    % get the channel at which max value is found
    %         [max_, maxORI, maxSF_log] = getIndMax(tempCovE, ORI_bound, plotFlag);
    maxORI=nan;maxSF_log=nan;
    max_allT(itrial, :) = [maxORI, maxSF_log];
    
    %     quickPlot_tempCovE, waitforbuttonpress
    % decide the format of IV
    switch IVType
        case 1, IV(itrial) = sum(tempCovE(:));
        case 2, IV(itrial) = max_;
        case 3, IV(itrial) = max_/(sum(tempCovE(:))); if sum(tempCovE(:)) == 0, error('ALERT: The suppressive drive is 0.'); end% SS stands for the estimated semi-saturation factor
    end
end % end of itrial

