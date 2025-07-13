%---------------------------------------------------------------
% Demo Script: Simulating Underestimation of Response Consistency
% Created by Shutian Xue on October 22, 2024
%
% This simulation illustrates why response consistency is often 
% underestimated in the trial-wise noisy observer model.
%---------------------------------------------------------------

% Setup
clear all; clc; close all;

% Add function paths
addpath(genpath('fxn_exp'))            % Experimental functions
addpath(genpath('fxn_analysis_RC_v2')) % Analysis functions

% Number of orientation channels
nORI = 29;

% Function to compute spatial-domain SD based on SF
fxn_getSigma_SPdomain = @(SF) 3 * sqrt(2*log(2)) / (2 * pi * SF); 

%--------------------------------------------%
% Load experimental settings
SX_RC1_setting
%--------------------------------------------%

% Define noise sampling range (used in filter construction)
noise.SF_low_sampling  = noise.SF_low;
noise.SF_high_sampling = noise.SF_high;

% Stimulus and noise parameters
cst_ln_template = 1.0;   % Contrast of the "true" Gabor template
noise.SF_low    = 1;     % Lower bound for SF sampling
noise.SF_high   = 4;     % Upper bound for SF sampling
gaborSD         = 0.8;   % Gabor SD in degrees

% Stimulus structure
ppd       = 32;           % Pixels per degree
sz_dva    = 3;            % Stimulus size (deg)
sz_pix    = sz_dva * ppd; % Stimulus size (pixels)
signalORI = 90;           % Orientation of the signal
gaborSF   = 2;            % SF of signal Gabor (cycles per deg)

stim.psz        = sz_dva;
stim.aper_psz   = sz_pix;
stim.targetOri  = signalORI;
stim.phase      = 0;
stim.gaborSF    = gaborSF;
stim.gaborSD    = gaborSD;
stim.gabor_sz   = sz_dva;
stim.mask       = exp_CreateCircularApertureSin(stim);

% Noise structure
noise.ppd                = ppd;
noise.psz                = sz_pix;
noise.fix_contrast       = 1;
noise.ratio_gaborInTgt   = 0.5; % Gabor-to-noise ratio in target
noise.ratio_base         = 0.5; % Base contrast ratio for noise

%% Generate Gabor filter bank (for visualization)
[filter_sin, filter_cos] = SX_sim02_setFilters( ...
    stim, noise.filtersSF_all, filtersOri_all, ...
    fxn_getSigma_SPdomain, 1);

fprintf('\nFilter bank created: nORI = %d, nSF = %d\n', ...
    length(filtersOri_all), length(noise.filtersSF_all));

%% Generate TRUE template (energy profile of the signal Gabor)

% Recompute filters (no plotting this time)
[filter_sin, filter_cos] = SX_sim02_setFilters( ...
    stim, noise.filtersSF_all, filtersOri_all, ...
    fxn_getSigma_SPdomain, 0);

close all;

% Choose SD for the input Gabor (normalized or fixed)
stim.gaborSD = fxn_getSigma_SPdomain(stim.gaborSF);
% stim.gaborSD = 4;  % Optionally override with fixed SD

% Create Gabor template patch
template_true_gabor = exp_CreateGabor(stim, cst_ln_template);

% Compute energy profile across filters
template_true = SX_RC4_Energy_parfor(stim.mask, ...
    {template_true_gabor}, filter_sin, filter_cos);
template_true = squeeze(template_true); % Remove singleton dim

% Normalize template values for visualization
new_max = 0.4;
new_min = -0.05;
template_true = (template_true - min(template_true(:))) * ...
    (new_max - new_min) / ...
    (max(template_true(:)) - min(template_true(:))) + new_min;

% Marginals
margORI_true = mean(template_true, 2);
margSF_true  = mean(template_true, 1);

% Plot: 2D Template + ORI and SF marginals
figure('Position', [0 0 2000 500])

% --- 2D Template ---
subplot(1, 3, 1); hold on;
imagesc(template_true);
[max_ORI, max_SF] = find(template_true == max(template_true(:)));
plot(max_SF, max_ORI, 'b*', 'MarkerSize', 10);
xline(15, 'r-');
yline(15, 'r-');
colorbar; clim([new_min, new_max]);
xlabel('SF channel #');
ylabel('ORI channel #');

% --- ORI Marginal ---
subplot(1, 3, 2); hold on;
plot(axis_tuning{1}, margORI_true, 'k-');
xticks(axisTicks_tuning{1});
yline(0, 'k--');
xline(0, 'r-');
xlabel('Orientation (º)');
ylim([new_min, new_max]);

% --- SF Marginal ---
subplot(1, 3, 3); hold on;
plot(axis_tuning{2}, margSF_true, 'k-');
xticks(axisTicks_tuning{2});
xticklabels(round(2.^axisTicks_tuning{2}, 2));
yline(0, 'k--');
xline(log2(2), 'r-');  % Reference line at 2 cpd
xline(axis_tuning{2}(margSF_true == max(margSF_true)), 'b-');
xlabel('SF (cpd)');
ylim([new_min, new_max]);

% --- Final Formatting ---
sgtitle(sprintf('True Template: Energy Profile of 2 cpd Gabor (Gabor SD = %.2f°)', stim.gaborSD));
set(findall(gcf, '-property', 'fontsize'), 'fontsize', 25);
set(findall(gcf, '-property', 'linewidth'), 'linewidth', 2);

% Optionally save figure
% saveas(gcf, sprintf('Fig/Sim_NOM_TrialWise/true_template_energy_profile_GaborSD%.2f.jpg', stim.gaborSD));