
% when use just one noise patch, rotate the patch to create different
% display

% if it's titration
if exp_mode == 0
    if noise.noiseVersion == 1 % one patch, but rotated
        [stimPatch, RotInd] = rotatePatch(noise.noisePatch{1}, size(gabor,1), lastRotInd);
        run.lastRotInd = RotInd;
    else % not rotated
        randn2 = randperm(noise.numNoisePatch);
        stimPatch = noise.noisePatch{randn2(1)}; % use the prepared noise patch
    end
    % if it's not titration
else, stimPatch = stim_noise; % generate a new noise patch in each trial
end