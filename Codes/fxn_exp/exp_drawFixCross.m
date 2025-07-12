function exp_drawFixCross(screen, stim)
centerX = screen.centerX;
centerY = screen.centerY;

fixSize = stim.fixSz;

Screen('DrawLine', screen.wPtr ,screen.black, centerX-fixSize, centerY, centerX+fixSize,centerY,  stim.PH_wid); 
Screen('DrawLine', screen.wPtr ,screen.black, centerX, centerY-fixSize, centerX, centerY+fixSize, stim.PH_wid);
