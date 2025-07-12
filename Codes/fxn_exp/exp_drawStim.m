
function run = exp_drawStim(exp_mode, params, run, record)

design = params.design;
stim = params.stim;
noise = params.noise;

limit0to1 = @(x) min(max(x,0),1);

if exp_mode ~= 1, run.contrast = run.stair.xCurrent; end

% decide to prs 5 or 1 stim (target)
if stim.stim5, nCuedLoc = design.nCuedLoc; else, nCuedLoc = 1;end

stimPatch_all = cell(nCuedLoc, 1);
noisePatch_all = stimPatch_all;

for iLoc = 1:design.nPrsLoc
    if ~stim.stim5, iLoc = run.iCuedLoc; end
    
    % decide contrast of gabor
    if exp_mode == 1, cst = record.thresh_exp(iLoc); else, cst = run.contrast; end
    
    % generate noisy gabor and noise for each loc
    [gabor, phase_tgt] = exp_CreateGabor(stim, cst, 0);
    filtered_noise = exp_CreateFilteredNoise(noise);
    stim_tgt = limit0to1(noise.ratio_base + gabor*noise.ratio_gaborInTgt + filtered_noise).*stim.mask + noise.ratio_base*(1-stim.mask);
    stim_noise = limit0to1(noise.ratio_base + filtered_noise).*stim.mask + noise.ratio_base*(1 - stim.mask);
    
    % decide the texture of stimulus
    if iLoc == run.iCuedLoc % the target
        if run.tgtPrs, stimPatch = stim_tgt; run.phase_tgt = phase_tgt;
        else, stimPatch = stim_noise; run.phase_tgt = nan; %   rotateStim_fxn % not in use
        end
    else % distractors (50% chance being a signal+noise, 50% chance being a pure noise)
        %         if randi(2,1) == 1, stimPatch = stim_tgt;
        %        else, stimPatch = stim_noise;
        %         end
        stimPatch = stim_noise;
    end
    
    stimPatch_all{iLoc} = stimPatch;
    noisePatch_all{iLoc} = filtered_noise;
end

run.stimPatch_all = stimPatch_all;
run.noisePatch_all = noisePatch_all;

