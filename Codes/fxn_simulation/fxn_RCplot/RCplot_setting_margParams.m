if nLoc==2
    % Params - ORI
    if paramMode==1, ticks{1, 1} = 0:.08:.32; else, ticks{1, 1} = 0:.08:.32; end % gain/peak amp
    ticks{1, 2} =  0:15:60; % width
    ticks{1, 3} = -.06:.04:.1; % baselin
    
    % Params - SF
    ticks{2, 1} = 0:.5:2; % peak SF (on log2 scale)
    if paramMode==1, ticks{2, 2} = 0:.06:.24; else, ticks{2, 2} = 0:.04:.16; end % gain vs. peak amp
    if flag_plotOctave, ticks{2, 3} = .4:.6:2.8; else, ticks{2, 3} = 1:.5:3.5; end % octave vs. cpd
    ticks{2, 4} = -.1:.05:.05; % baseline
    ticks{2, 5} = -.04:.04:.08; % truncation
    %     if flag_plotOctave, ydiffErr{2} = [1.08, .1, 2.1, .015, .055]; else, ydiffErr{2} = [1.08, .1, 1, .015, .055]; end
    
else % nLoc=4
    % Params - ORI
    ticks{1, 1} = 0:.05:.2; % gain/peak amp
    ticks{1, 2} =  5:10:35; % bandwidth
    ticks{1, 3} = -.05:.025:.05; % baseline
    % Params - SF
    ticks{2, 1} = 0:.5:1.5; % peak SF (on log2 scale)
    ticks{2, 2} = 0:.02:.1; % gain
    if flag_plotOctave, ticks{2, 3} = .5:.5:2.5; else, ticks{2, 3} = -.5:.5:1.5; end % width (on log2 scale)
    ticks{2, 4} = -.09:.03:.03; % baseline
    ticks{2, 5} = -.06:.04:.1; % truncation
end

