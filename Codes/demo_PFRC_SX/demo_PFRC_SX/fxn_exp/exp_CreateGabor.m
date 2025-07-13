function [gabor, phase] = exp_CreateGabor(stim, cst, plotQ)

% to generate Gabor with a circular window
if nargin==2, plotQ=0; end
[x1,y1] = meshgrid(1:stim.aper_psz, 1:stim.aper_psz); % size is the number of pixels
x2 = x1-mean(x1(:));
y2 = y1-mean(y1(:));
X = Scale(x2)*stim.aper_sz-stim.aper_sz/2; % Scale: perform an affine scaling to put data in range [0-1].
Y = Scale(y2)*stim.aper_sz-stim.aper_sz/2; % gabor_sz is sz*ppd

ori = stim.targetOri / 180 * pi; % all stim are always horizontally oriented
x_rotate = X*cos(ori) + Y*sin(ori); % The rotated coordinates
 
if ~isfield(stim, 'phase'), phase =  rand * 2 *pi; else, phase = stim.phase; end
sigma = stim.gaborSD;
carrier = cos(x_rotate * stim.gaborSF * 2 * pi + phase); % gaborSF_ppd is 0.0625
% modulator = exp(-((X / sigma).^2+(Y / sigma).^2)); % arch
modulator = exp(-(X.^2+Y.^2) / sigma^2); % new and correct
modulator = exp(-(X.^2+Y.^2) / (2*sigma^2)); % I missed a 2* in the deominator!!
gabor = carrier .* modulator * cst;

%% plot
if plotQ ==1
    flatSpreadVisRadius = (stim.gabor_psz(1).*stim.flatSpread)/2;
    figure('position',[0 600 1000 200]), clf, hold on,
    
    subplot(1,4,1)
    mesh(modulator)
    title('modulator')
    xlabel('X Pixel Loc')
    ylabel('Y Pixel Loc')
    zlabel('Contrast')
    
    subplot(1,4,2)
    mesh(carrier)
    title('carrier')
    xlabel('X Pixel Loc')
    ylabel('Y Pixel Loc')
    zlabel('Contrast')
    
    subplot(1,4,3)
    mesh(gabor)
    title('gabor')
    xlabel('X Pixel Loc')
    ylabel('Y Pixel Loc')
    zlabel('Brightness - ? Gray')
    hold off
    
    subplot(1,4,4)
    imshow(gabor)
    
    sgtitle(sprintf('SF = %.1f spread DVA = %d', stim.gaborSF*32,flatSpreadVisRadius))
end

    function output = Scale(input)
        % output = Scale(input)
        % Perform an affine scaling to put data in range [0-1].
        minval = min(input(:));
        maxval = max(input(:));
        output = (input - minval) ./ (maxval-minval);
    end
end