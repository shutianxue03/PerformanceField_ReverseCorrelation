
function [DV_allT, max_allT] = fxn_getDV_v3(e3D, template_input, convolveType, DV_allTType, flag_permT, ORI_bound)

flag_plotMax=0;

% Normalize template (Unit L2-norm)
template_input = template_input / norm(template_input(:));

if ndims(e3D)==3
    [nTrials, nORI, nSF] = size(e3D);
    DV_allT = nan(nTrials, 1);
    max_allT = nan(nTrials, 2);
else
    nTrials = 1;
    [nORI, nSF] = size(e3D);
    DV_allT = nan;
    max_allT = nan(1,2);
end

% Add Gaussian noise to the template, controlled by sigma_template.
% This is done by multiplying a Gaussian matrix (sigma = 1) with
% the energy profiles, then adding it to the integration between
% the true template and the energy profile.
% e3D_noisy = e3D .* randn(size(e3D));

% CANNOT permutate pixels of template here
% Because to have the same DV_allT for paired trials, the template cannot be randomized
% independently across trials (at least not for the pair of trials)

if flag_permT % Permute pixels in the template
    indRand = randperm(nORI*nSF);
    template_input_v = template_input(:);
    template_rand = template_input_v(indRand);
    template = reshape(template_rand, nORI, nSF);
else
    template = template_input;
end

for iTrial = 1:nTrials
    % This step is so important...
    if ndims(e3D)==3
        e2D = squeeze(e3D(iTrial, :, :));
    else
        e2D = e3D;
    end

    switch convolveType
        case 1, tempCovE = template.*e2D; % cross-correlation
        case 2 , tempCovE = conv2(template, e2D, 'same'); % convolution
    end

    % get the channel at which max value is found
    % [max_, maxORI, maxSF_log] = getIndMax(tempCovE, ORI_bound, flag_plotMax); % ORI_bound should be a vector of index of ori channels [ori_lb, ori_ub]
    % maxORI=nan;maxSF_log=nan;
    % max_allT(iTrial, :) = [maxORI, maxSF_log];

    %     quickPlot_tempCovE, waitforbuttonpress
    % decide the format of DV_allT
    switch DV_allTType
        case 1, DV_allT(iTrial) = sum(tempCovE(:));
        case 2, DV_allT(iTrial) = max_;
        case 3, DV_allT(iTrial) = max_/(sum(tempCovE(:))); if sum(tempCovE(:)) == 0, error('ALERT: The suppressDV_allTe drDV_allTe is 0.'); end% SS stands for the estimated semi-saturation factor
    end
end % end of itrial

