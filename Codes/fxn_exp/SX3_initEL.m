function EL = SX3_initEL(subjName, window, screen)

eyeDataDir = 'eyedata';
eyeFile = 'xx';
if length(eyeFile) > 8, error('EL file name CANNOT contain > 8 digits !!!'),end

% !!! file name cannot contain '_' !!!!
% file name can only contain 8 digits !!!

%% Initialize eye tracker
[EL,exitFlag] = rd_eyeLink('eyestart', window, {eyeFile, screen});
if exitFlag, return, end

%% Calibrate eye tracker
% [~, exitFlag] = rd_eyeLink('calibrate', window, EL);
% if exitFlag, return, end

EL.eyeDataDir  = eyeDataDir;
EL.eyeFile = eyeFile;

