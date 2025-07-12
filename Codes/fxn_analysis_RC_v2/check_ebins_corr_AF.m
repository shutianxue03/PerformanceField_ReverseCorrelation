%% check the correlation between energy bins
% each bin is the energy component of a specific ORI and SF channel of all
% trials

% 1. ExoEndo_DATASETS
cd('/Volumes/purplab/EXPERIMENTS/1_Current_Experiments/Antonio/Projects/ExoEndoRC_HPC_DATASETS/endo')
load('AS.mat')
load('exptParams.mat')
stiPar = expt.stiPar;
scr = expt.scr;
sz = stiPar.apersiz;
psz = stiPar.gaborpsiz;
sf = stiPar.gaborf;

% second we set up the filters
freStep          = logspace(log10(1.25),log10(3.25),15); nSF = length(freStep);
oriStep          = linspace(-80,80,19); nORI = length(oriStep);

xaxis_tuning = {oriStep, log2(freStep)};
xticks_tuning = {[-80, -44, 0, 44, 80], log2(freStep([1,4,8,12,15]))};
xticks_contour_tuning = {[1, 5, 10, 15, 19], [1,4,8,12,15]};

cueLocation = trialsMat(:,5); % trialsMat is the matrix that contains all of the observer's data
for ii = 1:length(cueLocation)
    if cueLocation(ii) == 1 || cueLocation(ii) == 3
        targettype(ii,:) = trialsMat(ii,3);
    elseif cueLocation(ii) == 2 || cueLocation(ii) == 4
        targettype(ii,:) = trialsMat(ii,4);
    end
end

%% get patch
icond = 1; namesCond = {'valid', 'neutral', 'invalid'};
iPRS = 0; text_PRSABS = ['ABS - ',  namesCond{icond}];
icst = 1; cst_unik = unique(trialsMat(:, 9)); 

% iPRS = 1; text_PRSABS = ['PRS - ',  namesCond{icond}];
patch_all = patch(targettype == iPRS & trialsMat(:, 2) == icond & trialsMat(:, 9) == cst_unik(icst));

% get energy
[e3D_AF, ~] = getEnergy(freStep, oriStep, patch_all, stiPar, scr);
fprintf('Energy calculated.\n')

%% reorganize AF's data (the format for check_corr_ebins is ntrials x nORI x nSF)
e3D_AF = dX;
[nSF, nORI, ntrials] = size(e3D_AF);
e3D_AF_flipped = nan(ntrials, nORI, nSF);
for itrial = 1:ntrials
    e3D_AF_flipped(itrial, :, :) = squeeze(e3D_AF(:, :, itrial))';
end

%%
function gabor = mkGabor(scr, gaborsiz, gaborstd, gaborf, phase, contrast, ori)

visiblesize=angle2pix(scr, [gaborsiz gaborsiz]);
[x,y]=meshgrid(1:visiblesize, 1:visiblesize);
x = x-mean(x(:));
y = y-mean(y(:));
x = Scale(x)*gaborsiz-gaborsiz/2;
y = Scale(y)*gaborsiz-gaborsiz/2;

ori = ori/180*pi;
nx = x*cos(ori) + y*sin(ori);
ny = x*cos(ori) - y*sin(ori);

carrier      =cos(nx*gaborf*2*pi+phase*2*pi);
modulator    =exp(-((x/gaborstd).^2)-((y/gaborstd).^2));
gabor        = carrier.*modulator*contrast;
end

%%
function [patch, full_objsiz] = mkMask(scr, objsiz, sinsiz, sinpower, imgsiz)

if nargin < 2 || isempty(objsiz)
    error('Not enough input arguments.');
end
if nargin < 3 || isempty(sinsiz)
    sinsiz = .5;
end
if nargin < 4 || isempty(sinpower)
    sinpower = 5;
end
if nargin < 5 || isempty(imgsiz)
    imgsiz = objsiz;
end
if objsiz > imgsiz
    error('Object size should be smaller than image size')
end

pimgsiz = angle2pix(scr, imgsiz);
psinsiz  = angle2pix(scr, sinsiz);

[x,y]=meshgrid(1:pimgsiz, 1:pimgsiz);
x = x-mean(x(:));
y = y-mean(y(:));
x = Scale(x)*imgsiz-imgsiz/2;
y = Scale(y)*imgsiz-imgsiz/2;

[~,r] = cart2pol(x,y);
patch = double(r <= objsiz/2);

sinFilter = sin(linspace(0,pi,psinsiz)).^sinpower;
sinFilter = sinFilter'*sinFilter;
sinFilter = sinFilter/sum(sinFilter(:));

patch = conv2(patch, sinFilter,'same');

Idx_1 = abs(patch-1) <= 10^-2;
Idx_r = r(Idx_1);
full_objsiz = 2*max(Idx_r(:));

end

%%
function pix = angle2pix(display,ang)

%Calculate pixel size
pixSize = display.width/display.resolution(1);   %cm/pix

sz = 2*display.dist*tan(pi*ang/(2*180));  %cm

pix = round(sz/pixSize);   %pix
end

%%
function [FR, PhaseValue] = getEnergy(freStep, OriStep, patch, stiPar, scr)
nFre = length(freStep);
nOri = length(OriStep);
nImg = length(patch);
FR = nan(nFre, nOri, nImg);
PhaseValue     = nan(nFre, nOri, nImg);
lumibg = patch{1}(1,1);

for i = 1:nFre
    fre = freStep(i);
    for j = 1:nOri
        ori = OriStep(j);
        
        filter1 = mkGabor(scr, stiPar.apersiz, stiPar.gaborstd, fre, 0, 1, ori);
        filter2 = mkGabor(scr, stiPar.apersiz, stiPar.gaborstd, fre, pi/2, 1, ori);
        
        mask  = mkMask(scr, stiPar.gaborsiz, stiPar.sinsiz, [], stiPar.apersiz);
        
        filter1 = filter1.*mask;
        filter2 = filter2.*mask;
        filter1 = filter1/sum(filter1(:).^2);
        filter2 = filter2/sum(filter2(:).^2);
        
        for k = 1:nImg
            tempImg = patch{k}-lumibg;
            sinImg = tempImg(:)'*filter1(:);
            cosImg = tempImg(:)'*filter2(:);
            FR(i,j,k)             = sqrt((sinImg)^2 + (cosImg)^2);
            PhaseValue(i,j,k)     = atan2(cosImg, sinImg);
        end
    end
end
end