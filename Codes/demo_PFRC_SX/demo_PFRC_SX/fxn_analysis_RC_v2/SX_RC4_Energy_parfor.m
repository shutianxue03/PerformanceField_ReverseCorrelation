function [e2D_allT, phase2D_allT] = SX_RC4_Energy_parfor(mask, patch_allT, filter_sin, filter_cos)
% SX_RC4_Energy_parfor: Computes energy and phase across trials and channels

% This function calculates the energy and phase for each orientation (ORI) and
% spatial frequency (SF) channel of each trial. 

% INPUTS:
%   mask            - the Gaussian mask 
%   patch_allT      - Cell array containing stimulus patches for all trials
%   filter_sin      - Cell array of sinusoidal filters (sin component), indexed by ORI and SF
%   filter_cos      - Cell array of sinusoidal filters (cos component), indexed by ORI and SF

% OUTPUTS:
%   e2D_allT        - 3D array storing energy values for each trial, ORI, and SF
%   phase2D_allT    - 3D array storing phase values for each trial, ORI, and SF

% Determine the number of orientations and spatial frequencies
[nORI, nSF] = size(filter_sin);
plot4parts = 0;  
nAllTrials = length(patch_allT);

% Initialize output arrays for energy and phase values across trials
e2D_allT = nan(nAllTrials, nORI, nSF);
phase2D_allT = e2D_allT;

% loop over orientations
for iORI = 1:nORI
    % fprintf('ORI #%d/%d...', iORI, nORI)
    
    % loop over spatial frequencies 
    parfor iSF = 1:nSF
        % Retrieve the sin and cos filters for the current orientation and SF
        filter_sin_ = filter_sin{iORI, iSF};
        filter_cos_ = filter_cos{iORI, iSF};
        
        % Loop over all trials to calculate energy and phase
        for itrial = 1:nAllTrials
            % Compute energy and phase for the current trial, orientation, and SF
            [energy, phase] = SX_sim04_computeEnergy(mask, patch_allT{itrial}, filter_sin_, filter_cos_, plot4parts);
            
            % Store results in the corresponding position in the output arrays
            e2D_allT(itrial, iORI, iSF) = energy;
            phase2D_allT(itrial, iORI, iSF) = phase;
        end
    end % End iSF loop
    
    % fprintf('DONE\n')
end % End iORI loop


%
% function [e2D_allT, phase2D_allT] = SX_RC4_Energy_parfor(mask, patch_allT, filter_sin, filter_cos)
%
% [nORI, nSF] = size(filter_sin);
% plot4parts = 0;     % 1 = plot the stim and two filters
% nAllTrials = length(patch_allT);
%
% % empty containers
% e2D_allT = nan(nAllTrials, nORI, nSF);
% phase2D_allT = e2D_allT;
%
% % loop to calculate energy of each channel
% for iORI = 1:nORI
% %     fprintf('ORI #%d/%d...', iORI, nORI)
%
%     parfor iSF = 1:nSF
%         filter_sin_ = filter_sin{iORI, iSF};
%         filter_cos_ = filter_cos{iORI, iSF};
%         for itrial = 1:nAllTrials
%
%             [a,b] = SX_sim04_computeEnergy(mask, patch_allT{itrial}, filter_sin_, filter_cos_, plot4parts);
%
%             e2D_allT(itrial, iORI, iSF) = a;
%             phase2D_allT(itrial, iORI, iSF) = b;
%         end
%
%     end % iSF
% %     fprintf('DONE\n')
% end % iORI
