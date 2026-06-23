%% Create three audience-demonstration videos
% Requires Image Processing Toolbox only for imresize; otherwise pure MATLAB.
% Outputs:
%   Figures/InteractiveDemoVideos/trial_fovea_*.mp4
%   Figures/InteractiveDemoVideos/trial_HM_*.mp4
%   Figures/InteractiveDemoVideos/trial_UVM_*.mp4
clear; close all; clc;
rng(20260626);

scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDir);
dirOutput = fullfile(projectRoot, 'Figures', 'InteractiveDemoVideos');
if ~exist(dirOutput, 'dir')
    mkdir(dirOutput);
end

%% Display and video settings
screenPx = [960, 1280];             % [height, width], matching the experiment
frameRate = 50;                     % Exact multiples for 400, 100, and 500 ms
background = 0.5;                   % Mid-gray
pxPerDeg = 40;                      % Adjustable for the defense display
eccDeg = 6;
patchDiameterDeg = 3;
fixationArmDeg = 0.3;
placeholderHalfWidthDeg = 1.5;
placeholderArmDeg = 0.45;
itiDuration = 0.000;
fixationDuration = 0.600;
stimulusDuration = 0.300; % increase to make the Gabor more visible for demonstration purposes
isiDuration = 0.300;
responseDuration = 1.200;

%% Stimulus settings from the paper
signalSF = 2;                       % cpd
signalSigmaSp = 6;                  % Larger values make the Gabor envelope broader
signalOrientation = 90;             % Horizontal grating
signalRMS = [.045];     % Titrated in the study; adjust for projector
noiseRMS = 0.10;
noiseLowSF = 1;
noiseHighSF = 4;
displayDistractors = 'off';         % 'on' or 'off'

%% Convert visual-angle parameters to pixels
eccPx = round(eccDeg * pxPerDeg);
patchDiameterPx = round(patchDiameterDeg * pxPerDeg);
fixationArmPx = round(fixationArmDeg * pxPerDeg);
placeholderHalfWidthPx = round(placeholderHalfWidthDeg * pxPerDeg);
placeholderArmPx = round(placeholderArmDeg * pxPerDeg);
if mod(patchDiameterPx, 2) == 0
    patchDiameterPx = patchDiameterPx + 1;
end
screenCenter = round([screenPx(2), screenPx(1)] / 2); % [x, y]
locations.fovea = screenCenter;
locations.rightHM = screenCenter + [eccPx, 0];
locations.leftHM = screenCenter + [-eccPx, 0];
locations.UVM = screenCenter + [0, -eccPx];
locations.LVM = screenCenter + [0, eccPx];
locationNames = {'fovea', 'rightHM', 'leftHM', 'UVM', 'LVM'};
locationXY = cell2mat(cellfun(@(name) locations.(name), locationNames, 'UniformOutput', false)');

%% Generate fixed distractors
%------------------------------------%
[~, aperture] = makeNoisyGaborPatch(patchDiameterPx, pxPerDeg, signalSF, signalSigmaSp, signalOrientation, signalRMS(1), noiseRMS, noiseLowSF, noiseHighSF);
%------------------------------------%
distractorPatches = cell(1, 5);
for iLoc = 1:5
    %----------------------------------------------%
    distractorPatches{iLoc} = makeFilteredNoisePatch(patchDiameterPx, pxPerDeg, noiseRMS, noiseLowSF, noiseHighSF, aperture);
    %----------------------------------------------%
end

%% Create videos for each signal RMS
conditions = {
    'fovea',   'fovea.mp4'
    'rightHM', 'HM.mp4'
    'UVM',     'UVM.mp4'
};
for iSignalRMS = 1:numel(signalRMS)
    thisSignalRMS = signalRMS(iSignalRMS);
    %--------------------------------------------%
    [targetPatch, aperture] = makeNoisyGaborPatch(patchDiameterPx, pxPerDeg, signalSF, signalSigmaSp, signalOrientation, thisSignalRMS, noiseRMS, noiseLowSF, noiseHighSF);
    %--------------------------------------------%
    rmsLabel = sprintf('%.1f', thisSignalRMS*100);
    for iCond = 1:size(conditions, 1)
        targetLocation = conditions{iCond, 1};
        [~, baseName, extension] = fileparts(conditions{iCond, 2});
        nameFile_video = fullfile(dirOutput, sprintf('%s_%s%s', baseName, rmsLabel, extension));
        itiFrame = background * ones(screenPx);
        itiFrame = drawPlaceholders(itiFrame, locationXY, locations.(targetLocation), placeholderHalfWidthPx, placeholderArmPx, fixationArmPx, false);
        fixationFrame = itiFrame;
        %----------------------------%
        fixationFrame = drawPlaceholders(fixationFrame, locationXY, locations.(targetLocation), placeholderHalfWidthPx, placeholderArmPx, fixationArmPx, true);
        %----------------------------%
        stimulusFrame = fixationFrame;
        distractorIndex = 1;
        for iLoc = 1:numel(locationNames)
            thisName = locationNames{iLoc};
            thisXY = locations.(thisName);
            if strcmp(thisName, targetLocation)
                patch = targetPatch;
                stimulusFrame = pastePatch(stimulusFrame, patch, aperture, thisXY);
            else
                if shouldDisplayDistractors(displayDistractors)
                    patch = distractorPatches{distractorIndex};
                    stimulusFrame = pastePatch(stimulusFrame, patch, aperture, thisXY);
                end
                distractorIndex = distractorIndex + 1;
            end
        end % iLoc
        isiFrame = itiFrame;
        isiFrame = drawPlaceholders(isiFrame, locationXY, locations.(targetLocation), placeholderHalfWidthPx, placeholderArmPx, fixationArmPx, false);
        responseFrame = isiFrame;
        %----------------------------%
        responseFrame = drawTextSimple(responseFrame, 'Was there a Gabor in the indicated patch?', round(screenPx(1) * 0.83));
        %----------------------------%
        %----------------------------%
        writeTrialVideo(nameFile_video, itiFrame, fixationFrame, stimulusFrame, isiFrame, responseFrame, frameRate, itiDuration, fixationDuration, stimulusDuration, isiDuration, responseDuration);
        %----------------------------%
        fprintf('Created %s\n', nameFile_video);
    end % iCond
end % iSignalRMS

%% Optional verification
fprintf('For each signal RMS, the same targetPatch matrix was reused across target locations.\n');

%% Local functions
function [patch, aperture] = makeNoisyGaborPatch(nPx, pxPerDeg, sfCpd, signalSigmaSp, orientationDeg, signalRMS, noiseRMS, lowSF, highSF)
    [xPx, yPx] = meshgrid(-(nPx-1)/2:(nPx-1)/2);
    xDeg = xPx / pxPerDeg;
    yDeg = yPx / pxPerDeg;
    theta = deg2rad(orientationDeg);
    xRot = xDeg * cos(theta) + yDeg * sin(theta);
    sigmaSp = signalSigmaSp * sqrt(2 * log(2)) / (2 * pi * sfCpd);
    phase = 2 * pi * rand;
    gabor = cos(2 * pi * sfCpd * xRot + phase) .* exp(-(xDeg.^2 + yDeg.^2) / (2 * sigmaSp^2));
    gabor = gabor - mean(gabor(:));
    gabor = gabor / std(gabor(:)) * signalRMS;
    aperture = makeSmoothCircularAperture(nPx, 0.45 * nPx, 0.05 * nPx);
    noise = makeBandpassNoise(nPx, pxPerDeg, lowSF, highSF, noiseRMS);
    contrastImage = gabor + noise;
    patch = 0.5 + contrastImage;
    patch = patch .* aperture + 0.5 * (1 - aperture);
    patch = min(max(patch, 0), 1);
end
function patch = makeFilteredNoisePatch(nPx, pxPerDeg, noiseRMS, lowSF, highSF, aperture)
    noise = makeBandpassNoise(nPx, pxPerDeg, lowSF, highSF, noiseRMS);
    patch = 0.5 + noise;
    patch = patch .* aperture + 0.5 * (1 - aperture);
    patch = min(max(patch, 0), 1);
end
function doDisplay = shouldDisplayDistractors(displayDistractors)
    doDisplay = strcmpi(displayDistractors, 'on');
end
function noise = makeBandpassNoise(nPx, pxPerDeg, lowSF, highSF, desiredRMS)
    whiteNoise = randn(nPx);
    frequencyAxis = ((0:nPx-1) - floor(nPx/2)) / nPx * pxPerDeg;
    [fx, fy] = meshgrid(frequencyAxis, frequencyAxis);
    radialSF = sqrt(fx.^2 + fy.^2);
    bandpassMask = radialSF >= lowSF & radialSF <= highSF;
    noiseSpectrum = fftshift(fft2(whiteNoise)) .* bandpassMask;
    noise = real(ifft2(ifftshift(noiseSpectrum)));
    noise = noise - mean(noise(:));
    noise = noise / std(noise(:)) * desiredRMS;
end
function aperture = makeSmoothCircularAperture(nPx, radiusPx, edgePx)
    [x, y] = meshgrid(-(nPx-1)/2:(nPx-1)/2);
    radius = sqrt(x.^2 + y.^2);
    aperture = 0.5 - 0.5 * tanh((radius - radiusPx) / edgePx);
end
function frame = pastePatch(frame, patch, aperture, centerXY)
    nPx = size(patch, 1);
    halfSize = floor(nPx / 2);
    xRange = centerXY(1)-halfSize:centerXY(1)+halfSize;
    yRange = centerXY(2)-halfSize:centerXY(2)+halfSize;
    validX = xRange >= 1 & xRange <= size(frame, 2);
    validY = yRange >= 1 & yRange <= size(frame, 1);
    xRange = xRange(validX);
    yRange = yRange(validY);
    patchX = find(validX);
    patchY = find(validY);
    localPatch = patch(patchY, patchX);
    localAperture = aperture(patchY, patchX);
    frame(yRange, xRange) = frame(yRange, xRange) .* (1 - localAperture) + localPatch .* localAperture;
end
function frame = drawPlaceholders(frame, locationXY, targetXY, halfWidthPx, armPx, fixationArmPx, showFixation)
    for i = 1:size(locationXY, 1)
        isTarget = isequal(locationXY(i, :), targetXY);
        if isTarget
            value = 1;
            thickness = 3;
        else
            value = 0;
            thickness = 2;
        end
        frame = drawCornerPlaceholder(frame, locationXY(i, :), halfWidthPx, armPx, value, thickness);
    end
    if showFixation
        frame = drawLine(frame, [size(frame, 2)/2-fixationArmPx, size(frame, 1)/2], [size(frame, 2)/2+fixationArmPx, size(frame, 1)/2], 0, 3);
        frame = drawLine(frame, [size(frame, 2)/2, size(frame, 1)/2-fixationArmPx], [size(frame, 2)/2, size(frame, 1)/2+fixationArmPx], 0, 3);
    end
end
function frame = drawCornerPlaceholder(frame, centerXY, halfWidthPx, armPx, value, thickness)
    halfWidth = round(halfWidthPx);
    xL = centerXY(1) - halfWidth;
    xR = centerXY(1) + halfWidth;
    yT = centerXY(2) - halfWidth;
    yB = centerXY(2) + halfWidth;
    frame = drawLine(frame, [xL, yT], [xL+armPx, yT], value, thickness);
    frame = drawLine(frame, [xL, yT], [xL, yT+armPx], value, thickness);
    frame = drawLine(frame, [xR, yT], [xR-armPx, yT], value, thickness);
    frame = drawLine(frame, [xR, yT], [xR, yT+armPx], value, thickness);
    frame = drawLine(frame, [xL, yB], [xL+armPx, yB], value, thickness);
    frame = drawLine(frame, [xL, yB], [xL, yB-armPx], value, thickness);
    frame = drawLine(frame, [xR, yB], [xR-armPx, yB], value, thickness);
    frame = drawLine(frame, [xR, yB], [xR, yB-armPx], value, thickness);
end
function frame = drawLine(frame, startXY, endXY, value, thickness)
    nSamples = ceil(max(abs(endXY - startXY))) + 1;
    x = round(linspace(startXY(1), endXY(1), nSamples));
    y = round(linspace(startXY(2), endXY(2), nSamples));
    for t = -floor(thickness/2):floor(thickness/2)
        for s = -floor(thickness/2):floor(thickness/2)
            xx = min(max(x + s, 1), size(frame, 2));
            yy = min(max(y + t, 1), size(frame, 1));
            frame(sub2ind(size(frame), yy, xx)) = value;
        end
    end
end
function frame = drawTextSimple(frame, textString, yPosition)
    rgb = repmat(grayDoubleToUint8(frame), 1, 1, 3);
    if exist('insertText', 'file') == 2
        rgb = insertText(rgb, [size(frame, 2)/2, yPosition], textString, 'AnchorPoint', 'Center', 'FontSize', 30, 'TextColor', 'black', 'BoxOpacity', 0);
        frame = rgbUint8ToGrayDouble(rgb);
    end
end
function writeTrialVideo(filename, itiFrame, fixationFrame, stimulusFrame, isiFrame, responseFrame, frameRate, itiDuration, fixationDuration, stimulusDuration, isiDuration, responseDuration)
    writer = VideoWriter(filename, 'MPEG-4');
    writer.FrameRate = frameRate;
    writer.Quality = 100;
    open(writer);
    writeRepeatedFrame(writer, itiFrame, round(itiDuration * frameRate));
    writeRepeatedFrame(writer, fixationFrame, round(fixationDuration * frameRate));
    writeRepeatedFrame(writer, stimulusFrame, round(stimulusDuration * frameRate));
    writeRepeatedFrame(writer, isiFrame, round(isiDuration * frameRate));
    writeRepeatedFrame(writer, responseFrame, round(responseDuration * frameRate));
    close(writer);
end
function writeRepeatedFrame(writer, grayFrame, nFrames)
    rgbFrame = repmat(grayDoubleToUint8(grayFrame), 1, 1, 3);
    for i = 1:nFrames
        writeVideo(writer, rgbFrame);
    end
end
function grayUint8 = grayDoubleToUint8(grayDouble)
    grayDouble = min(max(grayDouble, 0), 1);
    grayUint8 = uint8(round(grayDouble * 255));
end
function grayDouble = rgbUint8ToGrayDouble(rgbUint8)
    rgbDouble = double(rgbUint8) / 255;
    grayDouble = 0.2989 * rgbDouble(:, :, 1) + 0.5870 * rgbDouble(:, :, 2) + 0.1140 * rgbDouble(:, :, 3);
end