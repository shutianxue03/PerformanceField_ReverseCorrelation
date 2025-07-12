function EL = SX3_initEL(window, screen)

eyeDataDir = 'eyedata';
eyeFile = 'xx';
if length(eyeFile) > 8, error('EL file name CANNOT contain > 8 digits !!!'),end
% *** file name MUST NOT contain '_' ***
% file name can only contain 8 digits !!!

%% Initialize eye tracker
[EL, exitFlag] = rd_eyeLink('eyestart', window, {eyeFile, screen});
if exitFlag, return, end

%% save
EL.eyeDataDir  = eyeDataDir;
EL.eyeFile = eyeFile;
