function exp_drawPH(cue, stim, screen,iCuedLoc)

if stim.stim5, nPH = size(stim.ecc_box,1); else, nPH = 1;end

for iPH = 1:nPH % each PH
    if ~stim.stim5, iPH = iCuedLoc;end
    % decide the color of PH
    if ~stim.stim5
        PH_color = stim.noncue_color;
    else
        if cue && (iPH == iCuedLoc),PH_color = stim.cue_color;
        else,PH_color = stim.noncue_color;end
    end
    
    % draw corner by corner
    for ii = 1:4
        corner_loc = stim.PH_loc_box(ii,:);
        center(1) = corner_loc(1)*stim.PH_sz + screen.centerX + stim.ecc_box(iPH,1);
        center(2) = corner_loc(2)*stim.PH_sz + screen.centerY + stim.ecc_box(iPH,2);
        lines(1:2,:) = [center(1), center(2), center(1), center(2)-corner_loc(2)*stim.PH_leng; ...
            center(1), center(2), center(1)-corner_loc(1)*stim.PH_leng, center(2)];
        
        for n = 1:2
            Screen('DrawLine', screen.wPtr ,...
                repmat(PH_color,1,3), lines(n,1), lines(n,2), lines(n,3), lines(n,4), stim.PH_wid);
        end
    end
end
