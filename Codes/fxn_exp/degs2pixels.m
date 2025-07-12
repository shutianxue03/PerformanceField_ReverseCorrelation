function pixels = degs2pixels(screen, degs)
% screenRes - the resolution of the monitor,can be 1D or 2D
% screenSize - the size of the monitor in cm, can be 1D or 2D
% distance - the viewing distance in cm.
% degs - the amount of deg that should be transformed to a number of pixels

pixSize = screen.width./screen.resolution; %calculates the size of a pixel in cm
degPerPix = atand(pixSize./screen.dist);
pixPerDeg = 1./degPerPix;
pixels = round(pixPerDeg.*degs);
