
disp('checking phase:')

%% try filters of diff phase

phase_all = 0:pi/4:pi;
nPhase = length(phase_all);
energy_byPhase_allLoc = cell(3,nLoc);

for itype = 1:ntypes
    fprintf('target-%s trials ...\n', typeNames{itype})
    switch itype, case 1, trialInd = 1:ntrials;  case 2, trialInd = ntrials+1:nAllTrials;  case 3, trialInd = 1:nAllTrials; end
    
    for iLoc = 1:nLoc
        kernel_OnePhase = nan(nPhase ,nfilters);
        for iphase = 1:nPhase
            phase = phase_all(iphase);
            energy_onePhase = nan(length(trialInd), nfilters);
            for itrial = trialInd
                gaborPhase = data_both(itrial, iLoc, 5); gaborPhase(isnan(gaborPhase)) = 0;
                patch = patch_both{itrial, iLoc};
                lumiBG = patch(1,1);
                for ifilter = 1:nfilters
                    stim.gaborSF_ppd = noise.filterSF_all_ppd(ifilter);
                    stim.phase = phase;% + gaborPhase;
                    filter_ = exp_CreateGabor(stim, 1, 0) .* stim.mask;
                    filter = filter_/sqrt(sum(filter_(:).^2));
                    patch =  patch - lumiBG;
                    if itrial > ntrials, itrial_ = itrial-ntrials; else, itrial_ = itrial; end
                    energy_onePhase(itrial_, ifilter) = sum(patch(:).* filter(:));
                end
            end
            
            % derive kernels
            if itype < 3 % prs and abs
                e_mean = mean(energy_onePhase); % mean of each column (i.e., of each SF channel)
                e_sd = std(energy_onePhase,[],1);    % std of each column
                energy_norm = energy_onePhase - repmat(e_mean, ntrials,1) ./ repmat(e_sd, ntrials,1);
            else % both
                energy_norm = normEnergy(energy_onePhase);
            end
            
            kernel = SX_sim07_RC(filterSF_all, energy_norm, data_both(trialInd, iLoc, 2));
            kernel_OnePhase(iphase, :) = kernel;
        end
        
        energy_byPhase_allLoc{itype, iLoc} = kernel_OnePhase;
        fprintf('  Loc %d done.\n', iLoc)
    end
end

save(fileName_phase, 'energy_byPhase_allLoc')

