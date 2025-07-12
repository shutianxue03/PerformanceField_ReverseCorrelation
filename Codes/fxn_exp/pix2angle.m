function ang = pix2angle(screen,pix)

% . Inputs:
% . display.dist (distance from screen (cm))
% . display.width (width of screen (cm))
% . display.resolution (number of pixels of display in horizontal direction)
%
%  Output: ang (visual angle)
%
%Warning: assumes isotropic (square) pixels

%Written 11/1/07 gmb zre

%Calculate pixel size
pixSize = screen.width/screen.resolution(1);   %cm/pix
sz = pix*pixSize;  %cm (duh)
ang = 360/pi*atan(sz/(2*screen(1).dist));
