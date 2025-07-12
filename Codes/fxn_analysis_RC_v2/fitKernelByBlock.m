


%% get energy
lumiBG = patch_both{1}(1,1);
energy = nan(nTrialsPerBlock, nfilters);
for itrial = 1:nTrialsPerBlock
    for ifilter = 1:nfilters
        [a,b] = SX_sim04_computeEnergy(patch_both{itrial} , lumiBG, SFfilter_sin{ifilter}, SFfilter_cos{ifilter}, 0, itrial <= ntrials, filterSF_all(ifilter));
        energy(itrial,ifilter) = a;
    end
end

clear patch_both

energy_norm = [];
for n = 1:2 % tgt-prs and tgt-abs
    if n==1, energy_ = energy(1:nTrialsPerBlock/2, :); else, energy_ = energy(nTrialsPerBlock/2+1:end, :);end
    e_mean = mean(energy_); % mean of each column (i.e., of each SF channel)
    e_sd = std(energy_,[],1);    % std of each column
    energy_norm = [energy_norm; (energy_ - repmat(e_mean, nTrialsPerBlock/2,1)) ./ repmat(e_sd, nTrialsPerBlock/2,1)];
end

%% compute kernel
[kernel, var_exp] = SX_sim07_RC(nfilters, energy_norm, resp_both_B');
kernel_allLoc_allB(:, iLoc_B, ib_byLoc(iLoc_B)) = kernel;




