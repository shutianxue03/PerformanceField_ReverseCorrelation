
function [IV, e_noisy, max_allT] = fxn_getIV_v4(iModelA, e3D, convolveType, IVType, template_true, ORI_bound)

[ntrials, nORI, nSF] = size(e3D);
IV = nan(ntrials, 1);
e_noisy = IV;
max_allT = nan(ntrials, 2);

% CANNOT permutate pixels of template 
% to have the same IV for paired trials, the template cannot be randomized
% independently across trials (at least not for the pair of trials)

if (iModelA == 2) || (iModelA == 4) % randomize template on each trial
    rng('shuffle')
    indRand = randperm(nORI*nSF);
    template_true_v = template_true(:);
    template_rand = template_true_v(indRand);
    template = reshape(template_rand, nORI, nSF);
else
    template = template_true;
end

parfor itrial = 1:ntrials
    e = squeeze(e3D(itrial, :, :));
    
    % Add Gaussian noise to the template, controlled by sigma_template.
    % This is done by multiplying a Gaussian matrix (sigma = 1) with
    % the energy profiles, then adding it to the integration between
    % the true template and the energy profile.
    e_noisy(itrial) = sum(e .* randn(size(e)), 'all');
    
    switch convolveType
        case 1, tempCovE = template.*e; % dot multiply
        case 2 , tempCovE = conv2(template, e, 'same'); % convolve
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

