function exp_screenRefresh(screen)
centeredRect = CenterRectOnPointd([0 0 screen.size(1) screen.size(2)], screen.centerX, screen.centerY);
Screen('FillRect', screen.wPtr, screen.gray, centeredRect);
Screen('Flip', screen.wPtr);