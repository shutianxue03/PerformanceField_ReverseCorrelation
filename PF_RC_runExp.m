
% written by Shutian Xue
% initiated on 01/17/2021
% updated on 02/28/2021, added eyetracking component
% updated on 06/24/2021, finished RC analysis part, and quick titration at the  start of each session

%% clean, clear and prepare
warning('off')
clc
close all
clear all
format compact
commandwindow; % force the cursor to go automatically to command window
addpath(genpath('fxn_exp/'))
addpath(genpath('fxn_analysis_titration/'))
addpath(genpath('fxn_analysis_RC/'))
addpath(genpath('Data/'))
addpath(genpath('palamedes/'))
addpath(genpath('staircasecode_MJ/'))

Screen('Preference', 'SkipSyncTests', 2);
exp_modes = {'Titration','Experiment','Practice', 'Quick titration'};
EL_modes = {'OFF','ON'};

%% enter the info
clc
subjName = string(input('         >>> Enter initial: '));%'SX';
exp_mode= input('         >>> Enter exp mode (0=Titration, 1=Experiment, 2=Practice; 3=quick titration; 4=simulation): ');%8;
EL_mode =  input('         >>> Eyetracker on? (1=ON, 0=OFF): '); % 0=NO eyetracking
demoFlag = 0; if exp_mode==2, demoFlag=input('         >>> Show demo? (1=YES, 0=NO): '); end
screen_mode = 0; % 0 = full screen; 1 = small screen

clc
fprintf('Subject: %s\nMode: %s\nEye-tracker: %s\n',subjName, exp_modes{exp_mode+1},EL_modes{EL_mode+1})

% try
%% initialize params
params = SX1_initParams(screen_mode, exp_mode);

%% initialize eyelink
window = params.screen.wPtr;
if EL_mode, EL = SX3_initEL(subjName, window, params.screen); end
EL.ON = EL_mode;

%% pull out threshold (if available)
[fileName, record] = SX4_checkFile(exp_mode, subjName, params);

%% display instruction
SX6_Intro(params)

%% run the exp
record = SX7_run_AllBlocks(exp_mode, fileName, record, EL, params, demoFlag);

%% end exp
SX12_endExp

%% POST-EXP
SX13_postEXP

% catch
%     if EL.ON, rd_eyeLink('eyestop', params.screen.wPtr, {EL.eyeFile, EL.eyeDataDir}); end
%     PsychPortAudio('Close', params.soundParams.pahandle);
%     Screen('CloseAll');
%     error('ALERT: something wrong in the middle !!');
% end

clear all

