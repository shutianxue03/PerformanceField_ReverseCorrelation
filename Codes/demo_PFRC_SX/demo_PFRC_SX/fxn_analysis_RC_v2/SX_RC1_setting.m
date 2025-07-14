
% if isempty(dir('Data_OOD')), mkdir('Data_OOD'), end

% Directory of the Data folder on the server
nameFolder_Data = '/Volumes/purplab/EXPERIMENTS/1_Current_Experiments/Shutian_server/PF_RC/Data';
% Include Data, Data_OOD, Data_MC, Data_NOM*

nameFolder_Figures = '/Volumes/purplab/EXPERIMENTS/1_Current_Experiments/Shutian_server/PF_RC/Figures';

%%
nORI = 29;
nBins=10;
binStrategy_all = {'equal', 'algorithm', 'manual'};
eyeD_all = [1,1,0,1,1,1,1,1,0,0,1,0]; % 1=right eye dominant; 0=left eye dominant
flag_block200 = 0; % 1; all observers are forced to have 200 blocks; 0=no
markers_allSubj = {'o', 's', 'd', '^', 'v', '<', '>', '+', 'p', 'h', 'x', 'o', 's', 'd', '^'}; 
nMarkersMax = 11;  % Maximum number of different markers to use in plotting

%% Exp params
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
noise.SF_low = 1;  % Lower limit of SF sampling
noise.SF_high = 4;  % Upper limit of SF sampling
noise.ppd = ppd;
noise.sz=sz_dva;
noise.psz = sz_pix;
noise.fix_contrast = 1;
noise.ratio_gaborInTgt = 0.5;  % Ratio of Gabor signal in target
noise.ratio_base = 0.5;  % Base ratio for noise

params.stim=stim;
params.noise=noise;
% save('Data_OOD/params', 'params')

%%
nbins_e = 5; % number of energy bins
nitp = 2;% number of multiples of interpolated pints
nLoc8 = 8;
nLoc5 = 5;
nLoc2 = 2;
nlines = 2;
CI_ratio = .68;

threshPerf = .7;
dprime_theo = norminv(threshPerf)-norminv(1-threshPerf); % corresponding to 75% accuracy

%% get ORI/SF filters
% ORI filter
switch nORI
    case 29, fOri = [5:5:50, 60:10:90]; % nORI=29
    case 19, fOri = [5:10:45, 60:10:90]; % nORI=19
end
% fOri = [2, 5:5:50, 60:15:90]; % nORI=29
filtersOri_all = [-flip(fOri), 0, fOri] + 90; % must +90 !! otherwise the ori filters are 90 deg phased off
nORI = length(filtersOri_all);

% SF filter
nSF = nORI; % could generate 0.5: 17, 29, 37
noise.SF_low = 1;              
noise.SF_low = 1.1869; 
noise.SF_high = 2/noise.SF_low*2; noise.filtersSF_all_log = linspace(log2(noise.SF_low), log2(noise.SF_high), nSF); noise.filtersSF_all = 2.^noise.filtersSF_all_log;
cut_ORI = 1:nORI; cut_SF = 1:nSF;
nORI = length(cut_ORI);
nSF = length(cut_SF);

filtersOri_all = filtersOri_all(cut_ORI);
nSF = length(cut_SF);
filtersSF_all_log = noise.filtersSF_all_log(cut_SF);
filtersSF_all = noise.filtersSF_all(cut_SF);

noise.SF_low_sampling = noise.SF_low;
noise.SF_high_sampling = noise.SF_high;

%% ticks and labels
axis_tuning{1} = filtersOri_all - 90;
axis_tuning{2} = filtersSF_all_log;
axisTicks_tuning = {-90:45:90, linspace(log2(noise.SF_low), log2(noise.SF_high), 5)}; % ticks (SF is on log scale)
axisTL_tuning = {axisTicks_tuning{1}, round(2.^axisTicks_tuning{2}, 2)}; % label (SF is on linear scale)
axisLim = {[-99, 99], [-.1, 2.1]};

limit0to1 = @(x) min(max(x, 0), 1);
nTrialsPerSess = 100;  % Number of trials per session

%% index
combInd = [1,8; 6,7; 5,3; 2,4; 1,6;1,7];
ncomb = size(combInd, 1);
plotInd = [1:ncomb; ncomb+1:2*ncomb];
plotOne = 1; %1=save tuning curve plots in ONE figure; 0=save in separate figures
ncomb4 = 4;
ncomb2 = 2;
ncomb6 = 6;
ncomb8 = 8;

%% Model params (for fitting tunig curves)
% names of the MODEL
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
    'Raised DoG'};% M13

% names of the model PARAMs
namesParams_all = {
    {'gain', 'width', 'baseline'}, ...                                                    % M1 gaussian, to fit ori
    {'peak SF', 'gain', 'tuning', 'baseline'}, ...                             % M2 log parabola, AF
    {'peak SF', 'gain', 'bandwidth', 'baseline', 'truncation'}, ...              % M3 truncated log parabola, AF
    {'peak SF', 'gain', 'bandwidth', 'baseline', 'power'}, ...                 % M4 raised gaussian, LiBarbotCarrasco_2006
    {'peak SF', 'gain', 'bandwidth'}, ...                                             % M5 double exponential, JigoCarrasco2020
    {'peak SF', 'gain', 'bandwidth'}, ...                                             % M6 skewed gaussian, as suggested by AF
    {'peak SF', 'gain', 'bandwidth', 'baseline', 'power', 'truncate'}, ... % M7 truncated raised gaussian
    {'gain1', 'gain2', 'sigma1', 'sigma_r', 'baseline'}, ...                          % M8 difference of gaussians, to fit ori
    {'peak SF', 'gain', 'bandwidth', 'baseline'}, ... % M9: gaussian for SF
    {'gain', 'kappa', 'baseline'}, ... % M10:
    {'peak ORI', 'gain', 'kappa', 'baseline'}, ... % M11: Gaussian with the mean free to vary
    {'peakSF1', 'gain1', 'width1','base2', 'peakSF2', 'gain2', 'width2', 'base2'}, ... % M12: double peak
    {'gain1', 'gain2', 'sigma1', 'sigma_r', 'power', 'baseline'}};                      % M13, raised DoG

if ~exist('flag_standEnergy', 'var'), flag_standEnergy = 1; end

% lb and ub that are consistent across models
ORI_gain_lb = 1e-3; ORI_gain_ub = .2;
ORI_width_lb = 1e-3; ORI_width_ub = 40;
ORI_base_lb = -.2; ORI_base_ub = ORI_gain_ub;

SF_peak_lb = 1e-3; SF_peak_ub = 4;
SF_gain_lb = 1e-3; SF_gain_ub = .2;
SF_width_lb = 1e-3; SF_width_ub = .4; % .2 is arbitrary
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
    [ORI_gain_ub, ORI_width_ub, ORI_base_ub], ... % M11, Gaussian
    [SF_peak_ub, SF_gain_ub, SF_width_ub, SF_base_ub, SF_peak_ub, SF_gain_ub, SF_width_ub, SF_base_ub], ... % M12 double peak
    [ORI_gain_ub, 1, ORI_width_ub, 1e2, 10, ORI_base_ub]}; % M8, DoG 'gain1', 'gain2', 'sigma1', sigma_ratio, 'baseline'},

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
    [ORI_gain_lb, 0, ORI_width_lb, ORI_width_lb, -10, ORI_base_lb]}; % M13, raised DoG
% ============================================
fitMode = 2; % 1 = SSE, 2 = MLE;


%% Names
namesMetrics = {'dprime', 'criterion', 'pC', 'pHit', 'pFA', 'pA', 'pA1', 'pA0', 'CS'}; nmetrics = length(namesMetrics);
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

% MC
namesMCmode = {'10-CV', 'LOOCV',   'InfoCriterion'}; nMCmode = length(namesMCmode);
namesMCmode_long = {'10-fold cross validation (deviance is nLL)', ...
    'Leave-One-Out CV (deviance is nLL)', ...
    'Information criterion (deviance is nLL)'};
namesYLabel = {'Deviance', 'Deviance', 'IC'};

namesParamsMode = {'estP', 'tunC'};

% model
namesModelA = {'Core', 'RandTemp', 'IO-Core'};
% nModelsA = length(namesModelA);

namesModelB = {'Constant', 'Induced', 'noNoise', 'ConstantNoLapse', 'InducedNoLapse'};
% nModelsB = length(namesModelB);

namesParams2 = {'Thresh', 'mu [PRS]', 'sigma [PRS]','mu [ABS]', 'sigma [ABS]'};
nparams2 = length(namesParams2);

namesParamsModel_all = {{'Lapse rate', 'Constant noise', 'Criterion'}, {'Lapse rate', 'Induced noise', 'Criterion'}, {'Lapse rate', 'Criterion'}, ...
    {'Constant noise', 'Criterion'}, {'induced noise', 'Criterion'}};
% namesParamsNOM = {'Induced noise', 'Constant noise', 'Threshold'};
% namesParamsModel = {'Induced noise', 'Constant noise', 'Threshold'};
namesNormalityTest = {'Raw', 'Exp', 'Sqrt', 'boxcox'}; % reciprocal is deleted
nNormTests = length(namesNormalityTest);

namesConvolveType = {'dot product', 'convolution'}; nConvolveType = length(namesConvolveType);
namesIVType = {'sum all channels', 'channel with max IV'}; nIVType = length(namesIVType);

% plot
namesFeature_axis = {'Orientation (deg)', 'Spatial frequency (cpd)'};
% names of the estParams/tuningC
clear namesTunC_unit_perF
namesTunC_unit_perF{1,1} = {'Gain (a.u.)', 'Sigma (deg)', 'Baseline (a.u.)'};
namesTunC_unit_perF{1,2} = {'Pref ORI (deg)', 'Peak amp. (a.u.)', 'Bandwith (deg)', 'Baseline (a.u.)'};
namesTunC_unit_perF{2,1} = {'Peak SF (cpd)', 'Gain (a.u.)', 'Sigma (cpd)', 'Baseline (a.u.)'};
namesTunC_unit_perF{2,2} = {'Peak SF (cpd)', 'Peak amp. (a.u.)', 'Bandwith (octave)', 'Baseline (a.u.)'};
namesTunC_unit_perF{3,1} = {'Peak SF (cpd)', 'Gain (a.u.)', 'Sigma (cpd)', 'Baseline (a.u.)', 'Truncation (a.u.)'};
namesTunC_unit_perF{3,2} = {'Peak SF (cpd)', 'Peak amp. (a.u.)', 'Bandwith (octave)', 'Baseline (a.u.)', 'Truncation (a.u.)'};
namesTunC_unit_perF{8,1} = {'Gain 1 (a.u.)', 'Gain 2 (a.u.)', 'Sigma 1 (deg)', 'ratio (a.u.)', 'Baseline (a.u.)'};
namesTunC_unit_perF{8,2} = {'Pref ORI (deg)', 'peak amp. (a.u.)', 'trough ori (deg)', 'trough mag. (a.u.)', 'bandwidth (deg)', 'baseline (a.u.)'};
namesTunC_unit_perF{12,2} = {'peak SF1 (cpd)', 'peak amp. 1 (a.u.)', 'bandwidth 1 (octave)', 'peak SF2 (cpd)', 'peak amp. 2 (a.u.)', 'bandwidth 2 (octave)'};
namesTunC_unit_perF{13,2} = {'peak amp. (a.u.)', 'trough ori (deg)', 'trough mag, (a.u.)', 'bandwidth (deg)', 'baseline (a.u.)'};

%%
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
