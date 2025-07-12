patch_tgt = cell(nblocks, nTrialsPerBlock);

for ib = 1:nblocks
    load(dataFileDir(ib).name, 'record') % only load the record
    iLoc_B = record.iCuedLoc(ib,1);       % the iLoc of this block
    for itrial = 1:nTrialsPerBlock, patch_tgt{ib, itrial} = record.run_allTrials{itrial}.stimPatch_all{iLoc_B}; end % for new data
end


%% reorder patch_eachSubj given locations
patch_eachSubj = cell(nAllTrials, nLoc);

for iLoc = 1:nLoc % put location in order (center, L, U, R, D)
    orderInd = iCuedLoc == locAll(iLoc);
    patch_eachSubj(:,iLoc) = reshape(patch_tgt(orderInd,:), nAllTrials, 1);
end

%% reorganize so that PRS-trials are before the ABS-trials
patch_both = cell(nAllTrials, nLoc);
for iLoc = 1:nLoc
    % reorder patches
    patch_both(:, iLoc) = [patch_eachSubj(boolean(ii), iLoc); patch_eachSubj(boolean(1-ii), iLoc)];
end


%% get CI
CI_mean = cell(1,nLoc);
CI_var = CI_mean;
CI_mean_norm = CI_mean;
CI_var_norm = CI_mean;

for iLoc = 1:nLoc, RC_getCI, end


