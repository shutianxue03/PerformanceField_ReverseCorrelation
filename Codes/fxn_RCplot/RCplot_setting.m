

% if ~ exist('nsubj', 'var'), nsubj=1; isubj=1;end
% 
% if ~ exist('plotAve', 'var')
%     if nsubj == 1, plotAve = 0;
%         if ~ exist('subjName', 'var'), subjName=input('What is the subjName: ');end
%     else, plotAve = input('Given nSubj>1, plot group data (0=No, 1=YES): ');
%     end
%     if plotAve == 0, nsubj=1;end
% end

%% load params
load('Data_OOD/params1.mat')

stim = params.stim;
sz = stim.aper_sz;
threshPerf = .7;
dprime_theo = norminv(threshPerf)-norminv(1-threshPerf); % corresponding to 75% accuracy

%%  to take the weighted average based on trial numbers
if nsubj>1, ntrialsProp = nAllTrials'/sum(nAllTrials); else, ntrialsProp = 1; end

%%
polarAngCircle = linspace(0, 2*pi, 1e3);
polarAngFull = (0:90:360)/180*pi;
polarIndFull = [4,3,2,5,4]; %
subplot_locs = [5,4,2,6,8]; % iplot ind in a 3x3 grid, 5=center, 4=LHM, 2=UVM, 6=RHM, 8=LVM
nx_itpl = 3; % number of interpolated points between two data points when drawing a smooth fitting curve

polarAng = polarAngFull;
polarInd = polarIndFull;
% nAllTrials = ntrials*2;

%% plotting settings
faceAlpha = 0.3; % for error band
colorsType = {'r', 'b', 'k'}; % gabor-prs, abs, both
marks_allSubj_full = {'o', 's', 'd', '^','v',  '<', '>','p', 'h', '+', 'x', '-'}; % for each subj
marks_allSubj = marks_allSubj_full(1:nsubj); 
assert(length(marks_allSubj) == nsubj)

%% color
colors_comb = [
    0,0,0; ...,     % center; black
    0, .75, 0; ..., % Left: light green
    1, 0, 0,; ...,   % upper: red
    0, .35, 0; ..., % right: dark green
    0, 0, 1; ...,    % lower: blue
    0, .5, 0; ...,   % HM: green
    .5, 0, 1; ...,    % VM: purple
    .5, .5, .5];      % peri: darker grey

colors5Loc = colors_comb(1:5, :);
