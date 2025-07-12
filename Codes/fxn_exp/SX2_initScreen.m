function screen = SX2_initScreen(screen_mode)

screen = struct('rectPix',{[]}, 'dist', {[]},'size', {[]}, 'resolution', {[]},'calib_filename', {'0001_titchener_130226.mat'});

%% define some params
if screen_mode, screen.resolution = [800 600]  ;
else, screen.resolution = [1280 960];    
end
screen.rectPix = [0 0 screen.resolution];

screen.dist = 57; % default: 57
screen.size = [40 30]; % default: 40 30
screen.width = screen.size(1);
screen.height = screen.size(2);
screens = Screen('Screens');
screen.num = max(screens);
screen.white = WhiteIndex(screen.num);
screen.black = BlackIndex(screen.num);
screen.gray = round((screen.white+screen.black)/2);
screen.inc = screen.white-screen.gray;
screen.ratio_base = .5;
screen.bgColor = screen.ratio_base*255;
screen.pixelSz = pix2angle(screen,1);

% open the screen
AssertOpenGL;
wPtr = Screen('OpenWindow', 0, screen.bgColor, screen.rectPix);
screen.centerX = screen.resolution(1)/2;
screen.centerY = screen.resolution(2)/2;

screen.ppd = degs2pixels(screen,1); % pixel per degree
screen.ppd = 32; %screen.ppd(1);
screen.rad = 1.5 * screen.ppd; % acceptable fixation radius, in pixel
screen.wPtr = wPtr;
% Screen('BlendFunction', wPtr, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
%load('gammatable_L1_09012021.mat', 'GammaTable')
load('GammaTable_L1_03022023.mat', 'GammaTable')
if ~screen_mode, Screen('LoadNormalizedGammaTable', wPtr, GammaTable); end
fprintf('If the real experiment is ongoing, make sure the gamma table is loaded.\n')
% Priority(MaxPriority(wPtr));
screen.fd = Screen('GetFlipInterval',wPtr);


