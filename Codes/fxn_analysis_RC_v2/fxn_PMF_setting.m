
%% params
namesCCC = {'stair', 'const', 'manual', 'nonS' 'all'}; nccc = length(namesCCC);
perfThresh_all = [70]; nPerf = length(perfThresh_all);
iPerf_plot = 1; % plot which threshold? index the vector above
PMF_models = {'Logistic', 'CumNorm', 'Gumbel',  'Weibull'}; nModels = length(PMF_models);
% PMF_models = {'Logistic', 'CumNorm'}; nModels = length(PMF_models);
colors_allM = {'r', 'g', 'c', 'b'}; 

%% PMF fitting
analysisModes = [1,1; 1,0; 0,1; 0,0]; % left col: whether to bin data; right column: whether ti filter data (see OOD_fitPMF)

fit.nBins = 20;
fit.nBoot = 1e2; % bootstrapping, using PAL_PFML_BootstrapParametric() or PAL_PFML_BootstrapNonParametric()
fit.curveX = 10.^linspace(-3, 0, 1e3); % linear cst
fit.options.MaxIter = 1e12;
fit.options.MaxFunEvals = 1e12;
% range
fit.alphaguess = linspace(log10(.2), log10(.9), 1e3);
fit.betaguess = 10.^[-3:.1:3];
fit.gammaguess = .45:.01:.55;
fit.lambdaguess = 0:.01:.1;

fit.paramsFree = [1 1 1 1]; % 1=free to vary, 0=fixed
fit.guessLimits = [.45 .55];
fit.lapseLimits = [0 .1];
fit.nParams = 4; 
fit.searchGrid = struct('alpha', fit.alphaguess,'beta', fit.betaguess,'gamma',fit.gammaguess, 'lambda',fit.lambdaguess);

%% names
namesLoc9 = {'Fovea', '4 deg-Left', '4 deg-UVM', '4 deg-Right', '4 deg-LVM', '8 deg-Left', '8 deg-UVM', '8 deg-Right', '8 deg-LVM'};
namesLocCollapseHM = {'Fovea', '4 deg-Left', '4 deg-UVM', '4 deg-Right', '4 deg-LVM', '8 deg-Left', '8 deg-UVM', '8 deg-Right', '8 deg-LVM', ...
    '4 deg-HM', '8 deg-HM'};
namesINE = {'Equivalent noise (c^2)', 'Efficiency (%)'};
namesINE_short = {'Neq', 'Eff'};
namesLocComb = {'Fovea', 'Left', 'UVM', 'Right', 'LVM', 'HM', 'VM', 'Peri'};
namesAsym = {'Ecc effect', 'HVA', 'VMA', 'L/R diff'};

%% colors
g1 = .5;
g2 = .65;
g3 = .9;

colors_single = [0,0,0; ...,     % 1. fovea; black
    0, g2, 0; ..., % 2. LHM4: light green
    1, 0, 0,; ...,   % 3. UVM4: red
    0, g1, 0; ..., % 4. RHM4: dark green
    0, 0, 1; ...     % 5. LVM4: blue
    0, g3, 0;...    % 6. LHM8: lighter green
    1, .5, .5;...     % 7. UVM8: light red
    0, g2, 0;...   % 8. RHM8: light green
    .5, .75, 1];    % 9. LVM8: light blue

colors_comb = [
    ones(1,3)*0; ..., % 1. fovea: black
    ones(1,3)*.5; ..., % 2. ecc 4: grey
    ones(1,3)*.75; ...,   % 3. ecc 8: lighter grey
    0, g1, 0; ..., % 4. HM4: dark green
    .75, 0, .75; ...,    % 5. VM4: purple
    0, g2, 0; ...,   % 6. HM8: light green
    .75, .5, .75];      % 7. VM8: light purple

%% polar plots
% polar_ang = [180, 90, 0, 270]/180*pi;

%% 
iplots9 = [13, 12, 8, 14, 18, 11, 3, 15, 23];
iplots5 = [5, 4, 2, 6, 8];

%% compare across loc
iCompPair = [1,8; 6,7; 5,3; 2,4];
ncomp = size(iCompPair, 1);
getAsym = @(a,b) (a-b)./(a+b);

markers_allSubj = {'o', 's', 'd', '^','v',  '<', '+','p', 'h', 'x', '>', '*', 'p', 'o', 's', 'd', '^','v',  '<', '+','p', 'h', 'x', '>', '*', 'p', }; % for each subj

