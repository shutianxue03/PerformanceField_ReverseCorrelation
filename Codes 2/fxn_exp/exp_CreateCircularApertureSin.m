function [patch, full_objsiz] = exp_CreateCircularApertureSin(stim)

if nargin < 3 || isempty(stim.sin_sz)
    stim.sinsiz = .5;
end

if nargin < 4 || isempty(stim.sinpower)
    stim.sinpower = 5;
end

if nargin < 5 || isempty(stim.aper_sz)
    stim.aper_sz = stim.gabor_sz;
end

if stim.gabor_sz > stim.aper_sz
    error('Object size (gabor) should be smaller than image size (aperture).')
end

[x,y] = meshgrid(1:stim.aper_psz, 1:stim.aper_psz);
x = x-mean(x(:));
y = y-mean(y(:));
x = Scale(x)*stim.aper_sz-stim.aper_sz/2;
y = Scale(y)*stim.aper_sz-stim.aper_sz/2;

[~,r] = cart2pol(x,y);
patch = double(r <= stim.gabor_sz/2);

sinFilter = sin(linspace(0,pi,stim.sin_psz)).^stim.sinpower;
sinFilter = sinFilter'*sinFilter;
sinFilter = sinFilter/sum(sinFilter(:));

patch = conv2(patch, sinFilter,'same');

Idx_1 = abs(patch-1) <= 10^-2;
Idx_r = r(Idx_1);
full_objsiz = 2*max(Idx_r(:));

