if nLoc==2
    % Params - ORI
    yticks_all{1, 1} = -.3:.3:.6;
    yticks_all{1, 2} =  -.6:.4:.6; % width
    yticks_all{1, 3} = -.06:.04:.1; % baseline
    
    % Params - SF
    yticks_all{2, 1} = 0:.5:2; % peak SF (on log2 scale)
    yticks_all{2, 2} = 0:.16:.48; 
    if flag_plotOctave, yticks_all{2, 3} = .5:.5:2.5; else, yticks_all{2, 3} = 1:.5:3.5; end % octave vs. cpd
    yticks_all{2, 4} = -.1:.05:.05; % baseline
    yticks_all{2, 5} = -.04:.04:.08; % truncation
end

xticks_ = 0:.05:.2;
xticks_ = -.02:.05:.18;
% xticks_ = -.02:.03:.1;
xlim_ = xticks_ ([1, end]);