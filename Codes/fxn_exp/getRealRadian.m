function real_radian = getRealRadian(EL,mode,screen)
if mode==1
    EL.loc = Eyelink('NewestFloatSample'); % makes copy of most recent float sample received; returns -1 if no new sample or error
    
    if isempty(EL.loc) || ~isfield(EL.loc,'gx') || ~isfield(EL.loc,'gy')
        Eyelink('Message', sprintf('No fixation location reported. \n'));
        return
    else
        evt = Eyelink('newestfloatsample');
        EL.locX = evt.gx(2);%EL.guess.gx(2);
        EL.locY  = evt.gy(2);%EL.guess.gy(2);
    end
else
    mouseLoc = get(0, 'PointerLocation');
    EL.locX = mouseLoc(1);EL.locY = mouseLoc(2);
end

EL.locX = abs(EL.locX) - screen.centerX;
EL.locY = abs(EL.locY) - screen.centerY;

real_radian = sqrt((EL.locX^2) + (EL.locY^2));