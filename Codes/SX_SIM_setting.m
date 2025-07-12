%--------------%
SX_RC1_setting
%--------------%

nameFolder_NOM0 = 'Data_NOM_trialWise';

% Define the noise sampling range based on pre-set values
load(sprintf('Data_OOD/IO/signalEnergy_%d_%d.mat', nORI, nSF))
margORI_true = mean(template_true, 2);
margSF_true = mean(template_true, 1);

lapseRate=0;
ni = 1e3;
iModelA_all = [1]; 
iModelA_fit_all = 1;%[1, 5];  % 1: use the derived template, 5: use the ideal template
iModelB_fit_all = [4,5];
iLocComb=1;
patchMode = 1; % 1=energy calculated from target-patches; 2=from noise patches

%%
% Names of the performance metrics being analyzed
namesMetrics = {'pC', 'pYES', 'pA'};
nMetrics = length(namesMetrics);
nModelsA = 2; % do not use length()!!
nModelsB = 5;% do not use length()!!
nLocComb8 = 8;% do not use length()!!clc
clc
nBins = 10;  % Number of bins for the analysis
iIC_plot = 3; % % Which information criterion to plot (1 = AIC, 2 = AICc, 3 = BIC)

% Function handles for calculating information criteria (IC)
getAIC_SSE = @(SSE, nParams, nData) nData * log(SSE/nData) + 2*nParams;
getAICc_SSE = @(SSE, nParams, nData) nData * log(SSE/nData) + 2*nParams + 2*nParams*(nParams+1)/(nData-nParams-1);
getBIC_SSE = @(SSE, nParams, nData) nData * log(SSE/nData) + nParams * log(nData);

%%
% Define the noise sampling range based on pre-set values
noise.SF_low_sampling = noise.SF_low;
noise.SF_high_sampling = noise.SF_high;

% Initialize variables for the condition
cst_ln_template = 1.0;  % Constant for linear template
noise.noiseCST = noiseCST;
noise.SF_low = 1;  % Lower limit of SF sampling
noise.SF_high = 4;  % Upper limit of SF sampling
gaborSD = 0.8;  % Standard deviation of Gabor patch (default: 0.8)

% Define parameters related to stimulus
ppd = 32;  % Pixels per degree
sz_dva = 3;  % Size of stimulus in degrees of visual angle
sz_pix = sz_dva * ppd;  % Convert size to pixels
signalORI = 90;  % Orientation of signal (degrees)
gaborSF = 2;  % Spatial frequency of Gabor patches (cycles per degree)
stim.gaborCST = gaborCST;
stim.psz = sz_dva;
stim.aper_psz = sz_pix;
stim.targetOri = signalORI;
stim.phase = 0;
stim.gaborSF = gaborSF;
stim.gaborSD = gaborSD;
stim.gabor_sz = sz_dva;
stim.mask = exp_CreateCircularApertureSin(stim);

% Define noise parameters
noise.ppd = ppd;
noise.psz = sz_pix;
noise.fix_contrast = 1;
noise.ratio_gaborInTgt = 0.5;  % Ratio of Gabor signal in target
noise.ratio_base = 0.5;  % Base ratio for noise

limit0to1 = @(x) min(max(x, 0), 1);
nTrialsPerSess = 100;  % Number of trials per session
