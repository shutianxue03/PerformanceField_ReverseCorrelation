%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Script name: SX_RC1_setting.m
% Script type: script
% Author: Shutian Xue

% Last updated: 11/09/2025

% Description:
%   This script creates the settings for the PF-RC ANALYSIS pipeline.
%   It defines experimental settings (folder paths), global parameters (stimulus and noise properties), model fitting settings, names (in strings), and plotting styles.
%   The script is intended to be run at the start of most analysis scripts to ensure consistency and reproducibility across all modeling

% Define number of ORI (and SF) channels
nORI = 19;
nSF = nORI;

%% Define environment
nameFolder_current = pwd;

% Resolve project root: strip everything from the first /Codes component onward.
% This works whether the caller's pwd is PF_RC, PF_RC/Codes, or PF_RC/Codes/subfolder.
idx_codes = regexp(nameFolder_current, [filesep 'Codes' '($|' regexptranslate('escape', filesep) ')'], 'once', 'ignorecase');
if ~isempty(idx_codes)
    nameFolder_server = nameFolder_current(1:idx_codes - 1);
else
    nameFolder_server = nameFolder_current;
end

% Infer environment name from the project root path.
nameFolder_server_lower = lower(nameFolder_server);
if contains(nameFolder_server_lower, 'scratch')
    str_envir = 'HPC';
elseif contains(nameFolder_server_lower, 'server')
    str_envir = 'Server';
else
    str_envir = 'Local';
end

% nameFolder_server = '/Volumes/server/Users/purplab/EXPERIMENTS/1_Current_Experiments/Shutian_server/PF_RC';

%% Define names of folders to load/save data
% Keep Data/Figures parallel to Codes under the PF_RC project root.
nameFolder_Data = sprintf('%s/Data_Human', nameFolder_server) ;
nameFolder_Data_OOD = sprintf('%s/Data_OOD_%d%d', nameFolder_Data, nORI, nSF);  % Folder to save data
nameFolder_Data_NOM_Trialwise = sprintf('%s/Data_NOM_Trialwise_%d%d', nameFolder_Data, nORI, nSF);  % Folder to save data

% Define names of folders to save figures (on the server)
nameFolder_Figures = sprintf('%s/Figures', nameFolder_server);

% Define names of folders to save outputs (on the server)
nameFolder_Outputs = sprintf('%s/Outputs', nameFolder_server);

%% Define the function to normalize Gabor SD
fxn_getSigma_SPdomain = @(SF) 3 * sqrt(2*log(2)) / (2 * pi * SF);
% fxn_getSigma_SPdomain = @(SF) SF;

%%  Global parameters
nBins = 6; % to bin DVs in OOD_xx_beforeEst
binStrategy_all = {'equal', 'algorithm', 'manual'};
eyeD_all = [1,1,0,1,1,1,1,1,0,0,1,0]; % 1=right eye dominant; 0=left eye dominant
flag_block200 = 0; % 1; all observers are forced to have 200 blocks; 0=no
markers_allSubj_full = {'o', 's', 'd', '^', 'v', '<', '>', '+', 'p', 'h', 'x', 'o', 's', 'd', '^',     'o', 's', 'd', '^', 'v', '<', '>', '+', 'p', 'h', 'x', 'o', 's', 'd', '^',   'o', 's', 'd', '^', 'v', '<', '>', '+', 'p', 'h', 'x', 'o', 's', 'd', '^'};
nMarkersMax = 11;  % Maximum number of different markers to use in plotting
nbins_e = 5; % number of energy bins
nitp = 2;% number of multiples of interpolated pints
nLoc8 = 8;
nLoc5 = 5;
nLoc2 = 2;
nlines = 2;
CI_ratio = .68;

threshPerf = .7;
norminv_local = @(p) sqrt(2) * erfinv(2*p - 1);   % inverse standard-normal CDF (no toolbox required)
dprime_theo = norminv_local(threshPerf)-norminv_local(1-threshPerf); % corresponding to 75% accuracy

%% Stimulus and noise parameters
ppd = 32;  % Pixels per degree
sz_dva = 3;  % Size of stimulus in degrees of visual angle
sz_pix = sz_dva * ppd;  % Convert size to pixels
cst_ln_template = 1.0;  % Constant for linear template

% stim.gaborCST = gaborCST;
stim.ppd=ppd;
% stim.psz = sz_dva;
stim.aper_sz = sz_dva;
stim.aper_psz = sz_pix;
stim.gabor_sz = sz_dva;
stim.gabor_psz = sz_pix;
stim.sin_sz = .8;
stim.sinpower = 5;
stim.sin_psz = sz_dva;
stim.targetOri = 90;
stim.phase = 0;
stim.gaborSF = 2;% Spatial frequency of Gabor patches (cycles per degree)
stim.gaborSD = 0.8;% Standard deviation of Gabor patch (default: 0.8)
stim.mask = exp_CreateCircularApertureSin(stim);

% Define noise parameters
noise.noiseCST = .2;
noise.SF_low = 1;  % Lower limit %               6 = estimate multiplicative noise (Nmul) and SDadd;of SF sampling
noise.SF_high = 4;  % Upper limit of SF sampling
noise.ppd = ppd;
noise.sz=sz_dva;
noise.psz = sz_pix;
noise.fix_contrast = 1;
noise.ratio_gaborInTgt = 0.5;  % Ratio of Gabor signal in target
noise.ratio_base = 0.5;  % Base ratio for noise

params.stim=stim;
params.noise=noise;

%% Settings for ORI/SF filter
fOri = linspace(0, 90, (nORI+1)/2); fOri = fOri(2:end);
filtersOri_all = round([-flip(fOri), 0, fOri] + 90); % must +90 !! otherwise the ori filters are 90 deg phased off
assert(nORI == length(filtersOri_all));

% SF filter
nSF = nORI; % could generate 0.5: 17, 29, 37
noise.SF_low = 1; % used to be 1.1869
noise.SF_high = 2/noise.SF_low*2;
noise.filtersSF_all_log = linspace(log2(noise.SF_low), log2(noise.SF_high), nSF);
noise.filtersSF_all = 2.^noise.filtersSF_all_log;

filtersSF_all_log = noise.filtersSF_all_log;
filtersSF_all = noise.filtersSF_all;

noise.SF_low_sampling = noise.SF_low;
noise.SF_high_sampling = noise.SF_high;

%% Ticks and labels
axis_tuning{1} = filtersOri_all - 90;
axis_tuning{2} = filtersSF_all_log;
axisTicks_tuning = {-90:45:90, linspace(log2(noise.SF_low), log2(noise.SF_high), 5)}; % ticks (SF is on log scale)
axisTL_tuning = {axisTicks_tuning{1}, round(2.^axisTicks_tuning{2}, 2)}; % label (SF is on linear scale)
axisLim = {[-99, 99], [.1719, 1.8281]};

limit0to1 = @(x) min(max(x, 0), 1);
nTrialsPerSess = 100;  % Number of trials per session

%% Index
combInd = [1,8; 6,7; 5,3; 2,4; 1,6;1,7];
ncomb = size(combInd, 1);
plotInd = [1:ncomb; ncomb+1:2*ncomb];
plotOne = 1; %1=save tuning curve plots in ONE figure; 0=save in separate figures
ncomb4 = 4;
ncomb2 = 2;
ncomb6 = 6;
ncomb8 = 8;

%% NOM fitting
% Parameter bounds and initial values
NOMp1_lb = 0;
NOMp2_lb = 0;
NOMp3_lb  = 0;

NOMp1_ub = 2;     % multi noise
NOMp2_ub = 4;   % (private) additive noise
NOMp3_ub = 3; % (shared) additive noise

% Midpoint initial guesses
NOMp1_0 = mean([NOMp1_lb,  NOMp1_ub]);
NOMp2_0 = mean([NOMp2_lb, NOMp2_ub]);
NOMp3_0   = mean([NOMp3_lb,   NOMp3_ub]);

% BADS options
options_bads = bads('defaults');
options_bads.Display  = 'none';
options_bads.Verbosity = 0;

% fmincon options
options_fmin = optimoptions('fmincon', 'MaxIterations', 1e4, 'Display', 'off');

% Noisy observer model names
namesModelA = {'RC', 'IO', 'RandTemp'};

namesModelB = {'FullModel', 'NoMulti', 'NoAdd', 'NoShared', 'MultiOnly', 'AddOnly', 'SharedOnly', 'JustCriterion'};
flag_regressType = 2;  % 1=Univariate; 2=Multi+smoothing
flag_incluCrit   = 1;  % 1=include criterion as a free parameter
C_contribution   = 0;  % contribution of criterion loss to the objective

switch flag_incluCrit
    case 0
        namesModelBparams = {...
            {'Multiplicative variability', 'Additive variability', 'Shared variability'}, ... %1
            {                                         'Additive variability', 'Shared variability'}, ... %2
            {'Multiplicative variability',                                'Shared variability'}, ... %3
            {'Multiplicative variability', 'Additive variability'                                       }, ... %4
            {'Multiplicative variability',                                                                           }, ... %5
            {                                         'Additive variability',                                      }, ... %6
            {                                                                        'Shared variability'}, ... %7
            {                                                                                                  }}; %8
    case 1
        namesModelBparams = {...
            {'Multiplicative variability', 'Additive variability', 'Shared variability', 'CriterionNOM'}, ... %1
            {                                         'Additive variability', 'Shared variability',  'CriterionNOM'}, ... %2
            {'Multiplicative variability',                                'Shared variability',    'CriterionNOM'}, ... %3
            {'Multiplicative variability', 'Additive variability'                                , 'CriterionNOM'}, ... %4
            {'Multiplicative variability',                                                                  'CriterionNOM'}, ... %5
            {                                         'Additive variability',                                 'CriterionNOM'}, ... %6
            {                                                                        'Shared variability',    'CriterionNOM'}, ... %7
            {                                                                                                  'CriterionNOM'}}; %8

        namesModelBparams_short = {...
            {'Nmul', 'Nadd', 'Nshared', 'criterion_DV'}, ... %1
            {            'Nadd', 'Nshared',  'criterion_DV'}, ... %2
            {'Nmul',             'Nshared',  'criterion_DV'}, ... %3
            {'Nmul', 'Nadd'                  , 'criterion_DV'}, ... %4
            {'Nmul',                               'criterion_DV'}, ... %5
            {             'Nadd',                  'criterion_DV'}, ... %6
            {                       'Nshared',    'criterion_DV'}, ... %7
            {                                    'criterion_DV'}}; %8
otherwise
        error('Invalid value for flag_incluCrit. Must be 0 or 1.');
end
namesConvolveType = {'dot product', 'convolution'}; nConvolveType = length(namesConvolveType);
namesDVType = {'sum all channels', 'channel with max DV'}; nDVType = length(namesDVType);

%% Tuning function params
% names of the tuning models
namesFamily_all = {
    'Gaussian', ... % M1
    'log parabola', ... % M2
    'truncated parabola', ... % M3
    'raised gaussian', ... % M4
    'double exponential', ... % M5
    'skewed gaussian', ... % M6
    'truncated raised gaussian', ... % M7
    'DoG', ... % M8
    'Gaussian (SF)', ...
    'von Mises', ... % M10
    'Gaussian (mean free)', ...  % M11
    'Double peak', ...% M12
    'Raised DoG', ... % M13
    'Asymmetric Gaussian (SF)'};% M14

% names of the model PARAMs
namesParams_all = {
    {'gain', 'width', 'baseline'}, ...                                                    % M1 gaussian, to fit ori
    {'peak SF', 'gain', 'tuning', 'baseline'}, ...                             % M2 log parabola, AF
    {'peak SF', 'gain', 'width', 'baseline', 'truncation'}, ...              % M3 truncated log parabola, AF
    {'peak SF', 'gain', 'width', 'baseline', 'power'}, ...                 % M4 raised gaussian, LiBarbotCarrasco_2006
    {'peak SF', 'gain', 'width'}, ...                                             % M5 double exponential, JigoCarrasco2020
    {'peak SF', 'gain', 'width'}, ...                                             % M6 skewed gaussian, as suggested by AF
    {'peak SF', 'gain', 'width', 'baseline', 'power', 'truncate'}, ... % M7 truncated raised gaussian
    {'gain1', 'gain2', 'sigma1', 'sigma_r', 'baseline'}, ...                          % M8 difference of gaussians, to fit ori
    {'peak SF', 'gain', 'width', 'baseline'}, ... % M9: gaussian for SF
    {'gain', 'kappa', 'baseline'}, ... % M10:
    {'peak ORI', 'gain', 'width', 'baseline'}, ... % M11: Gaussian with the mean free to vary
    {'peakSF1', 'gain1', 'width1','base2', 'peakSF2', 'gain2', 'width2', 'base2'}, ... % M12: double peak
    {'gain1', 'gain2', 'sigma1', 'sigma_r', 'power', 'baseline'}, ... % M13, raised DoG
    {'peak SF', 'gain', 'width left', 'width right', 'baseline'}}; % M14, asymmetric Gaussian (SF)

if ~exist('flag_standEnergy', 'var'), flag_standEnergy = 1; end

% lb and ub that are consistent across models
ORI_gain_lb = 1e-3; ORI_gain_ub = .5;
ORI_width_lb = 1e-3; ORI_width_ub = 60;
ORI_base_lb = -.2; ORI_base_ub = ORI_gain_ub;

SF_peak_lb = 1e-3; SF_peak_ub = 4;
SF_gain_lb = 1e-3; SF_gain_ub = .5;
SF_width_lb = 1e-3; SF_width_ub = 1; % 1 is arbitrary; but 0.4 is too low
SF_base_lb = -.1; SF_base_ub = SF_gain_ub;

ub_full_all = {
    [ORI_gain_ub, ORI_width_ub, ORI_base_ub], ... % M1. Gaussiam [gain, width, baseline]
    [SF_peak_ub, SF_gain_ub, SF_width_ub, SF_base_ub], ... % M2. Log parabola {'peak SF', 'gain', 'tuning', 'baseline'}
    [SF_peak_ub, SF_gain_ub, SF_width_ub, SF_base_ub, SF_gain_ub], ... % M3 log parabola with truncation {'peak SF', 'gain', 'tuning', 'baseline', truncation'}
    [2, .2, 2, 1, 10], ... % M4 raised gaussian
    [2, .2, 2], ... % M5 double exponential
    [3.5, .2, 4], ... % M6 skewed gaussian
    [3.5, .2, 4, 1, 10, 10],... % M7 raised gaussian with truncation
    [ORI_gain_ub, 1, ORI_width_ub, 1e2, ORI_base_ub], ... % M8, DoG 'gain1', 'gain2', 'sigma1', sigma_ratio, 'baseline'},
    [2, 1, 2, 1], ...
    [.5, 10, .2], ... % M10, von Mises, {gain, kappa, base}
    [20, ORI_gain_ub, ORI_width_ub, ORI_base_ub], ... % M11, Gaussian with free mean [peak, gain, width, baseline]
    [SF_peak_ub, SF_gain_ub, SF_width_ub, SF_base_ub, SF_peak_ub, SF_gain_ub, SF_width_ub, SF_base_ub], ... % M12 double peak
    [ORI_gain_ub, 1, ORI_width_ub, 1e2, 10, ORI_base_ub], ... % M13, raised DoG
    [SF_peak_ub, SF_gain_ub, SF_width_ub, SF_width_ub, SF_base_ub]}; % M14, asymmetric Gaussian (SF)

lb_full_all = {
    [ORI_gain_lb, ORI_width_lb, ORI_base_lb], ... % M1. Gaussiam [gain, width, baseline]
    [SF_peak_lb, SF_gain_lb, SF_width_lb, SF_base_lb], ... % M2. Log parabola {'peak SF', 'gain', 'tuning', 'baseline'}
    [SF_peak_lb, SF_gain_lb, SF_width_lb, SF_base_lb, SF_base_lb], ... % M3 log parabola with truncation {'peak SF', 'gain', 'tuning', 'baseline', truncation'}
    [2, .2, 2, 1, 10], ... % M4 raised gaussian
    [2, .2, 2], ... % M5 dolble exponential
    [3.5, .2, 4], ... % M6 skewed gaussian
    [3.5, .2, 4, 1, 10, 10],... % M7 raised gaussian with truncation
    [ORI_gain_lb, 0, ORI_width_lb, ORI_width_lb, ORI_base_lb], ... % M8, DoG 'sigma', 'sigma_d', 'gain', 'baseline'},
    [.0001, .0001, 1, -.2], ...
    [.0001, .0001, -.2],... % M10, von Mises, {gain, kappa, base}
    [-20,ORI_gain_lb, ORI_width_lb, ORI_base_lb], ... % M11
    [SF_peak_lb, SF_gain_lb, SF_width_lb, SF_base_lb, SF_peak_lb, SF_gain_lb, SF_width_lb, SF_base_lb], ... % M12
    [ORI_gain_lb, 0, ORI_width_lb, ORI_width_lb, -10, ORI_base_lb], ... % M13, raised DoG
    [SF_peak_lb, SF_gain_lb, SF_width_lb, SF_width_lb, SF_base_lb]}; % M14, asymmetric Gaussian (SF)
% ============================================
fitMode = 2; % 1 = SSE, 2 = MLE;

%% Names
% General names
% metrics = [dprime, criterion, [pC, pHit, pFA], nanmean(respC), pYES, mean(1./contrast), median(RT)];
% namesMetrics = {'dprime', 'criterion', 'pC', 'pHit', 'pFA', 'pA', 'pA1', 'pA0', 'CS'}; nmetrics = length(namesMetrics);
namesMetrics = {'dprime', 'criterion', 'pC', 'pHit', 'pFA', 'pA', 'pA1', 'pA0', 'pYES', 'CS', 'RT'}; nMetrics = length(namesMetrics);
namesMetricsLong = {'Dprime', 'criterion', 'Accuracy', 'pHit', 'pFA', 'pA', 'pA1', 'pA0', 'pYES', 'Contrast sensitivity', 'Resp. time'}; nMetrics = length(namesMetrics);
namesFeature = {'ORI', 'SF'}; nFeatures = length(namesFeature);
namesType = {'PRS', 'ABS', 'BOTH'}; nTypes = length(namesType);
namesLoc2D = {'Fovea', 'Left', 'Upper', 'Right', 'Lower'};
namesLocComb = {'Fovea', 'LHM', 'UVM', 'RHM', 'LVM', 'HM', 'VM', 'Peri.'};
namesSDT = {'d''', 'criterion', 'RT (sec)'};
namesEnergy = {'WHOLE', 'NOISE', 'NOISE-masked'};
namesNorm = {'Unnorm', 'Norm'};
namesRespType = {'Hit rate', 'FA rate', 'pC'};
namesTitles1 = {'Fovea', 'Horizontal', 'Lower', 'Left'}; % before 'vs.', color-coded
namesTitles2 = {'Perifovea', 'Vertical', 'Upper', 'Right'}; % after 'vs.', color-coded
namesTitles3 = {'', 'Meridian', 'Vertical Meridian', 'HM'}; % in black
% namesCI = {'mean', 'mean norm', 'var', 'var norm'};
namesIC = {'AIC', 'AICc', 'BIC'}; nICs = length(namesIC);
publishOptions = struct('format','pdf','outputDir','publishedPDFs/', 'showCode', 0);
% namesDataset = {'TrainingSet', 'FullSet'};
% namesDataset = {'TempSet', 'FullSet'}; % template set, full set (all data)
namesDataset_full = {'FullSet', 'TmplSet', 'TrainSet', 'TestSet'}; nDatasets_full = length(namesDataset_full);
namesDataset = {'TmplSet', 'FullSet'}; nDatasets = length(namesDataset); % needs to matchOOD_xx_compDV ("for iDataset = 1:2")

% Model comparison names
namesMCmode = {'10-CV', 'LOOCV',   'InfoCriterion'}; nMCmode = length(namesMCmode);
namesMCmode_long = {'10-fold cross validation (deviance is nLL)', ...
    'Leave-One-Out CV (deviance is nLL)', ...
    'Information criterion (deviance is nLL)'};
namesYLabel = {'Deviance', 'Deviance', 'IC'};

namesParamsMode = {'estP', 'tunC'};

%% For plotting
namesFeature_axis = {'Orientation (º)', 'Spatial frequency (cpd)'};
namesFeature_axis_Tuning = {'Marginalized ORI weights (a.u.)', 'Marginalized SF weights (a.u.)'};
% names of the estParams/tuningC
clear namesTunC_unit_perF
namesTunC_unit_perF{1,1} = {'Gain (a.u.)', 'Sigma (deg)', 'baseline (a.u.)'};
% namesTunC_unit_perF{1,2} = {'Pref ORI (deg)', 'amplitude (a.u.)', 'width (deg)', 'baseline (a.u.)'};
namesTunC_unit_perF{1,2} = {'amplitude (a.u.)', 'width (º)', 'baseline (a.u.)'};
namesTunC_unit_perF{2,1} = {'peak SF (cpd)', 'Gain (a.u.)', 'Sigma (cpd)', 'baseline (a.u.)'};
namesTunC_unit_perF{2,2} = {'peak SF (cpd)', 'amplitude (a.u.)', 'width (octaves)', 'baseline (a.u.)'};
namesTunC_unit_perF{3,1} = {'peak SF (cpd)', 'Gain (a.u.)', 'Sigma (cpd)', 'baseline (a.u.)', 'Truncation (a.u.)'};
namesTunC_unit_perF{3,2} = {'peak SF (cpd)', 'amplitude (a.u.)', 'width (octaves)', 'baseline (a.u.)', 'Truncation (a.u.)'};
namesTunC_unit_perF{8,1} = {'Gain 1 (a.u.)', 'Gain 2 (a.u.)', 'Sigma 1 (deg)', 'ratio (a.u.)', 'baseline (a.u.)'};
namesTunC_unit_perF{8,2} = {'Pref ORI (deg)', 'amplitude (a.u.)', 'trough ori (deg)', 'trough mag. (a.u.)', 'width (deg)', 'baseline (a.u.)'};
namesTunC_unit_perF{10,1} = {'Gain (a.u.)', 'kappa (a.u.)', 'baseline (a.u.)'};
namesTunC_unit_perF{10,2} = {'amplitude (a.u.)', 'width (º)', 'baseline (a.u.)'};
namesTunC_unit_perF{11,1} = {'peak ORI (deg)', 'Gain (a.u.)', 'Sigma (deg)', 'baseline (a.u.)'};
namesTunC_unit_perF{11,2} = {'amplitude (a.u.)', 'width (º)', 'baseline (a.u.)'};
namesTunC_unit_perF{12,2} = {'peak SF1 (cpd)', 'amplitude 1 (a.u.)', 'width 1 (octaves)', 'peak SF2 (cpd)', 'amplitude 2 (a.u.)', 'width 2 (octaves)'};
namesTunC_unit_perF{13,2} = {'amplitude (a.u.)', 'trough ori (deg)', 'trough mag, (a.u.)', 'width (deg)', 'baseline (a.u.)'};
namesTunC_unit_perF{14,2} = {'peak SF (cpd)', 'amplitude (a.u.)', 'width (octaves)', 'baseline (a.u.)'};

namesTunC_noUnit{1,1} = {'Gain', 'Sigma', 'baseline'};
namesTunC_noUnit{1,2} = {'amplitude', 'width', 'baseline'};
namesTunC_noUnit{2,1} = {'peak SF', 'Gain', 'Sigma', 'baseline'};
namesTunC_noUnit{2,2} = {'peak SF', 'amplitude', 'width', 'baseline'};
namesTunC_noUnit{3,2} = {'peak SF', 'amplitude', 'width', 'baseline', 'Truncation'};
namesTunC_noUnit{8,2} = {'Pref ORI', 'amplitude', 'trough ori', 'trough mag.', 'width', 'baseline'};
namesTunC_noUnit{10,1} = {'Gain', 'kappa', 'baseline'};
namesTunC_noUnit{10,2} = {'amplitude', 'width', 'baseline'};
namesTunC_noUnit{11,1} = {'peak ORI', 'Gain', 'Sigma', 'baseline'};
namesTunC_noUnit{11,2} = {'amplitude', 'width', 'baseline'};
namesTunC_noUnit{12,2} = {'peak SF1', 'amplitude 1', 'width 1', 'peak SF2', 'amplitude 2', 'width 2'};
namesTunC_noUnit{13,2} = {'amplitude', 'trough ori', 'trough mag.', 'width', 'baseline'};
namesTunC_noUnit{14,2} = {'peak SF', 'amplitude', 'width', 'baseline'};


% Colors
colors_comb = [
    0,0,0; ...,     % center; black
    0, .75, 0; ..., % Left: light green
    0, 0, 1,; ...,   % upper: blue
    0, .35, 0; ..., % right: dark green
    1, 0, 0; ...,    % lower: red
    0, .5, 0; ...,   % HM: green
    .5, 0, 1; ...,    % VM: purple
    .5, .5, .5];      % peri: darker grey

colorsType = {'r', 'b', 'k'};
% colors2 = {'r', 'b'};
colors2 = {[1,0,0], [0,0,1]}; % when only present fovea and peri