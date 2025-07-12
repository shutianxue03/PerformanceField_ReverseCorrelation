%% sound
InitializePsychSound(1); % 1 for precise timing
% Open the audio device
pahandle = PsychPortAudio('Open');
device = PsychPortAudio('GetDevices');
params.fs = device.DefaultSampleRate;
params.pahandle =  pahandle;
params.freqCorrect = 1e3;
params.freqNeutral = 8e2;
params.freqWrong = 6e2;
params.amp = 1;
for ff = 6e2:2e2:1e3, makeBeep(ff, .2), end