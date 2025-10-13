function [energy, phaseValue] = SX_sim04_computeEnergy(mask, patch, filter_sin_raw, filter_cos_raw, plot4parts)
% Computes the energy and phase of a stimulus patch
% 
% This function calculates the energy and phase of an input stimulus patch by 
% applying pool of Gabor filters with quadrature phase (sin and cos)
% 
% INPUTS:
%   mask            - Gaussian mask applied to filters
%   patch           - 2D array of the stimulus patch
%   filter_sin_raw  - Raw sinusoidal filter component (sin)
%   filter_cos_raw  - Raw sinusoidal filter component (cos)

% OUTPUTS:
%   energy          - contrast energy 
%   phaseValue      - Phase angle (in radians)

% Exit function if patch is empty
if isempty(patch), return, end

% Define the background luminance (the first pixel of the patch)
lumiBG = patch(1,1);

% Apply mask to filters and normalize to unit energy
filter_sin_m = filter_sin_raw .* mask;
filter_cos_m = filter_cos_raw .* mask;
% filter_sin = filter_sin_m / sum(filter_sin_m(:).^2);   % Normalize sin filter
% filter_cos = filter_cos_m / sum(filter_cos_m(:).^2);   % Normalize cos filter
filter_sin = filter_sin_m / sqrt(sum(filter_sin_m(:).^2));   % Normalize sin filter by the length
filter_cos = filter_cos_m / sqrt(sum(filter_cos_m(:).^2));   % Normalize cos filter

% Remove background luminance from the patch
patch = patch - lumiBG;

% Calculate the response to the sin and cos filters
stim_sin = patch(:)' * filter_sin(:);  % Dot product with sin filter
stim_cos = patch(:)' * filter_cos(:);  % Dot product with cos filter

% Compute energy as the magnitude of the response vector
energy = sqrt(stim_sin^2 + stim_cos^2);

% Compute phase value using the arctangent of the sin and cos responses
phaseValue = atan2(stim_sin, stim_cos);

% Optional: Plot the stimulus and filter components if plot4parts is true
if plot4parts
    figure('Position', [0 300 400 400])
    subplot(2, 2, 1), imshow(patch), title('Stimulus Patch')
    subplot(2, 2, 2), imshow(filter_sin), title('Sin Filter (Normalized)')
    subplot(2, 2, 3), imshow(filter_sin.^2), title('Sin Filter Squared')
    subplot(2, 2, 4), imshow(filter_sin .* patch), title('Filter x Stimulus')
    suptitle(sprintf('Energy = %.4f\n', energy))
    pause(0.5)
end



% function [energy, phaseValue] = SX_sim04_computeEnergy(mask, patch, filter_sin_raw, filter_cos_raw, plot4parts)
% % FernandezLiCarrasco_2019
% if isempty(patch), return, end
% 
% lumiBG = patch(1,1);
% 
% filter_sin_m = filter_sin_raw.*mask;
% filter_cos_m = filter_cos_raw.*mask;
% filter_sin = filter_sin_m/sum(filter_sin_m(:).^2);   % normalize the filter
% filter_cos = filter_cos_m/sum(filter_cos_m(:).^2); % normalize the filter
% 
% patch =  patch - lumiBG;
% stim_sin = patch(:)' * filter_sin(:);
% stim_cos = patch(:)' * filter_cos(:);
% energy = sqrt(stim_sin^2 + stim_cos^2);
% phaseValue = atan2(stim_sin, stim_cos);
% 
% %%
% if plot4parts
%     figure('Position',[0 300 400 400])
%     subplot(2,2,1), imshow(patch),title('stim')
%     subplot(2,2,2), imshow(filter_sin),title('sin filter')
%     subplot(2,2,3), imshow(filter_sin.^2),title('sin filter^2')
%     subplot(2,2,4), imshow(filter_sin.*patch),title('filter x stim')
%     suptitle(sprintf('energy = %.4f\n',energy))
%     pause(.5)
% end
