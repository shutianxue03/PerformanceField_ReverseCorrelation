function makeBeep(freq, dur, amp)
global params
pahandle = params.pahandle;
fs  = params.fs;
if nargin==2, amp = params.amp; end

% create sound data
nBuffer = 100;
beep_  = amp*sin(2*pi* freq*(0:1/fs:dur));
beep_ = [linspace(0,beep_(1), nBuffer), beep_, linspace(beep_(end), 0, nBuffer)];

% Fill the audio playback buffer with the audio data 'wavedata':
PsychPortAudio('FillBuffer', pahandle, [beep_; beep_]);

% Start audio playback
PsychPortAudio('Start', pahandle);
PsychPortAudio('Stop', pahandle, 3);

