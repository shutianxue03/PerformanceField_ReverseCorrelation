
% this function organize important data
% outputs:
%       - data_eachSubj: nAllTrials x nLoc x 4 measures (itgtPrs, answer, correctness, RTs), location in order
%       - patch_eachSubj: nAllTrials x nLoc cell, each cell contains a patch; location in order
%       - data_both: rows reordered so that prs-trials are before abs-trials
%       - patch_both: rows reordered so that prs-trials are before abs-trials

%% load params
nSess = nblocks/nLoc;
nSess = 1; % mute this!!
nAllTrials = nSess * nTrialsPerBlock;
ntrials = nAllTrials/2;

%% extract patches (from each data file)
patch_tgt = cell(4, nTrialsPerBlock); % change 4 back to nblocks!
for ib = (nblocks-3):nblocks%% change back to 1:nblocks !!
    load(dataFileDir(ib).name, 'record') % only load the record
%     for itrial = 1:nTrialsPerBlock, patch_tgt{ib, itrial} = record.patch_allTrials{itrial}{record.iCuedLoc(ib, 1)}; end
    ind4 = mod(ib,4); if ind4==0, ind4 = 4;end % mute this!!
    for itrial = 1:nTrialsPerBlock, patch_tgt{ind4, itrial} = record.patch_allTrials{itrial}{record.iCuedLoc(ib, 1)}; end % mute this!!
end

% reorder given locations
patch_eachSubj = cell(nAllTrials, nLoc);
% iCuedLoc = record.iCuedLoc(1:nblocks, 1);
iCuedLoc = record.iCuedLoc(nblocks-3:nblocks, 1); % mute this!
for iLoc = 1:nLoc % put location in order (center, L, U, R, D)
    orderInd = iCuedLoc == locAll(iLoc);
    patch_eachSubj(:,iLoc) = reshape(patch_tgt(orderInd,:), nAllTrials, 1);
end

%% extract behavioral measurements (from the last data file)
nmeasures = 4; % itgtPrs, answer, correctness, RT
data_eachSubj = nan(nAllTrials, nLoc, nmeasures); % each layer is a set of data

for n = 1 : nmeasures
    switch n
        case 1, x_record = record.targetPrs;
        case 2, x_record = record.answer;
        case 3, x_record = record.correctness;
        case 4, x_record = record.RT;
    end
%     x_record = x_record(1:nblocks, 1:nTrialsPerBlock); % only take the non-nan part
    x_record = x_record(nblocks-3:nblocks, 1:nTrialsPerBlock); % mute this !!
    for iLoc = 1:nLoc % adjust the order according to the location
        orderInd = iCuedLoc == locAll(iLoc);
        data_eachSubj(:, iLoc, n) = reshape(x_record(orderInd, :), nAllTrials, 1);
    end
end

%% reorganize so that PRS-trials are before the ABS-trials
data_both = nan(nAllTrials, nLoc, nmeasures);
patch_both = cell(nAllTrials, nLoc);
for iLoc = 1:nLoc
    ii = data_eachSubj(:,iLoc, 1); % 1st layer is PrsInd
    for n = 1 : nmeasures
        data_both(:, iLoc, n) = [data_eachSubj(boolean(ii), iLoc, n); data_eachSubj(boolean(1-ii), iLoc, n)];
    end
    patch_both(:, iLoc) = [patch_eachSubj(boolean(ii),iLoc); patch_eachSubj(boolean(1-ii),iLoc)];
end

%% SDT
dprime_allLoc = nan(1, nLoc);
criterion_allLoc = dprime_allLoc;
RT_allLoc = dprime_allLoc;

for iLoc = 1:nLoc
    prs_both = data_both(:,iLoc, 1);
    ans_both = data_both(:,iLoc, 2);
    RT_both = data_both(:,iLoc, 4);
    
    nHit = sum(ans_both & prs_both);
    nFA = sum(ans_both & 1-prs_both);
    [dprime, criterion] = SX_sim06_SDT(nHit, nFA, ntrials);
    
    dprime_allLoc(iLoc) = dprime;
    criterion_allLoc(iLoc) = criterion;
    RT_allLoc(iLoc) = median(RT_both);
end