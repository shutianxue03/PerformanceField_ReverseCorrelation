
simMode = 1;
nLoc = 1;
iLoc = 1;

%% experiment params
% locName_all = {'Center','LHM','RHM','LVM', 'UVM'};
% nLoc = length(locName_all);
params = SX1_initParams(0,4);

noise = params.noise;
stim = params.stim;

screen.ratio_base = .5; % the BG luminance, slight brighter than the exp
screen.ppd = 32;
screen.resolution = [1280 960];    
screen.dist = 57; % default: 57
screen.size = [40 30]; % default: 40 30
screen.width = screen.size(1);
screen.height = screen.size(2);
noise.scr = screen;

%%



%% Names
locNames = {'Fovea', 'Left', 'Upper', 'Right', 'Lower'};
locAll = 1:length(locNames);
typeNames = {'prs', 'abs', 'both'}; 
ntypes = length(typeNames);
SDT_titles = {'d''', 'criterion', 'RT (sec)'};
publishOptions = struct('format','pdf','outputDir','publishedPDFs/', 'showCode', boolean(0));

% nLoc = params.design.nCuedLoc; assert(length(locNames) == nLoc)

%% Model params (for fitting tunig curves)
modelSF = 2; % the model used to fit SF kernels
modelOri = 1; % the model used to fit ori kernels
fitMode = 1; % 1 = SSE, 2 = MLE;
modelNames = {'Gaussian', 'log parabola', 'truncated parabola', 'raised gaussian', 'double exponential', 'skewed gaussian', 'truncated raised gaussian', 'DoG'};

kernelNames = {'SF', 'Ori'};

modelNameOri = modelNames{modelOri};
modelNameSF = modelNames{modelSF};

modelParamsNamesAll = {
    {'gain', 'width', 'baseline'}, ...                                                    % gaussian, to fit ori
    {'peak SF', 'gain', 'tuning', 'baseline'}, ...                             % log parabola, AF
    {'peak SF (cpd)', 'sigma', 'alpha', 'baseline', 'truncate'}, ...              % truncated parabol, AF
    {'peak SF (cpd)', 'sigma', 'alpha', 'baseline', 'power'}, ...                 % raised gaussian, LiBarbotCarrasco_2006
    {'peak SF (cpd)', 'sigma', 'alpha'}, ...                                             % double exponential, JigoCarrasco2020
    {'peak SF (cpd)', 'sigma', 'alpha'}, ...                                             % skewed gaussian, as suggested by AF
    {'peak SF (cpd)', 'sigma', 'alpha', 'baseline', 'power', 'truncate'}, ... % truncated raised gaussian
    {'sigma1', 'sigma2', 'gain1', 'gain2', 'baseline1', 'baseline2'}};                          % differnce of gaussians, to fit ori

modelParamsNamesOri = modelParamsNamesAll{modelOri};
modelParamsNamesSF = modelParamsNamesAll{modelSF};
nparamsOri = length(modelParamsNamesOri);
nparamsSF = length(modelParamsNamesSF);

%% get ORI/SF filters
nbins_e = 4; % number of energy bins

% ORI filter
fOri = [5, 7.5, 10, 12.5, 15, 20, 30, 40, 60, 80];
filtersOri_all = [-flip(fOri), 0, fOri] + 90;
nfiltersOri = length(filtersOri_all);

% SF filter
filtersSF_all = noise.filtersSF_all;
filtersSF_all_log = log2(filtersSF_all);
nfiltersSF = length(filtersSF_all);
% noise.filtersSF_all = filtersSF_all;





