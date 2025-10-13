
%% [Behav] measurements of each session
dprime_perSess_perLoc = nan(nSess, nLoc5);
criterion_perSess_perLoc = dprime_perSess_perLoc;
RT_perSess_perLoc = dprime_perSess_perLoc;
cst_perSess_perLoc = dprime_perSess_perLoc;
pA_perSess_perLoc = nan(ntypes, nSess, nLoc5);
pC_perSess_perLoc = pA_perSess_perLoc; % do not delete, for plotting pC vs. pA

for iLoc = 1:nLoc5
    dprime_perSess_perLoc(:, iLoc) = dprime_perBlk(iCuedLoc_allT==iLoc);
    criterion_perSess_perLoc(:, iLoc) = criterion_perBlk(iCuedLoc_allT==iLoc);
    RT_perSess_perLoc(:, iLoc) = RT_aveT_perBlk(iCuedLoc_allT==iLoc);
    cst_perSess_perLoc(:, iLoc) = cst_perBlk(iCuedLoc_allT == iLoc);
    pA_perSess_perLoc(:, :, iLoc) = pA_perBlk(:, iCuedLoc_allT==iLoc);
    pC_perSess_perLoc(:, :, iLoc) = pC_perBlk(:, iCuedLoc_allT==iLoc);
end

%% get the pass index of each trial (1=pass A, 2=pass B)
passInd = nan(nblocks, nTrialsPerBlock);

for iblk = 1:nblocks
    for ip = 1:npairs
        ind_samePass = find(pairInd(iblk, :) == ip);
        passInd(iblk, ind_samePass(1)) = 1;
        passInd(iblk, ind_samePass(2)) = 2;
    end
end

%% [Behav] reorder behavioral measurements given locations
nm7 = 7; % (1) itgtPrs, (2) answer, (3) correctness, (4) RT, (5) cst (6) ind of pass (1=A, 2=B) (7) pairInd within each pass
% data_eachSubj = nan(nAllTrials, nLoc, nmeasures); % each layer is a set of data
data_allT_perComb = cell(nLoc8, nm7);

for im = 1:nm7
    switch im
        case 1, x_record = tgtPrs_allT;
        case 2, x_record = answer_allT;
        case 3, x_record = correctness_allT;
        case 4, x_record = log(RT_allT);
        case 5, x_record = cst_allT;
        case 6, x_record = passInd;
        case 7, x_record = pairInd;
    end
    
    if size(x_record, 1) > nblocks % cst_allT has the correct cut shape
        x_record = x_record(ib_start:ib_end, :);
    end
    
    for iLoc = 1:nLoc8 % adjust the order according to the location
        if iLoc>5, switch iLoc, case 6, iLoc_all = [2,4]; case 7, iLoc_all = [3,5]; case 8, iLoc_all = 2:5; end
        else, iLoc_all = iLoc;
        end
        dd = cell(length(iLoc_all),1); % empty containers
        for ii = 1:length(iLoc_all)
            dd{ii} = reshape(x_record(iCuedLoc_allT == iLoc_all(ii), :), ntrials*2, 1);
        end
        data_allT_perComb{iLoc, im} = cell2mat(dd);
    end
end

%% [Patch] reorder patches given locations
if flagGetPatch
    patch_both_perComb_p = cell(nLoc8, 1);
    noise_both_perComb_p = patch_both_perComb_p;
    noiseMasked_both_perComb_p = patch_both_perComb_p;
    
    for ip = 1:3
        switch ip
            case 1, pp_allT = patch_tgt_perBlk;
            case 2, pp_allT = noise_tgt_perBlk;
            case 3, pp_allT = noiseMasked_tgt_perBlk;
        end
        
        for iLoc = 1:nLoc8 % put location in order (center, L, U, R, D)
            if iLoc>5, switch iLoc, case 6, iLoc_all = [2,4]; case 7, iLoc_all = [3,5]; case 8, iLoc_all = 2:5; end
            else, iLoc_all = iLoc;
            end
            
            pp = cell(length(iLoc_all)*ntrials, 1);
            for ii = 1:length(iLoc_all)
                pp(1+nAllTrials*(ii-1):nAllTrials*ii) = reshape(pp_allT(iCuedLoc_allT == iLoc_all(ii), :), nAllTrials, 1);
            end
            
            % reorder given itgt
            pp_reorder = cell(nAllTrials * length(iLoc_all),1);
            pp_reorder(1:length(iLoc_all)*ntrials) = pp(boolean(data_allT_perComb{iLoc, 1}));
            pp_reorder(length(iLoc_all)*ntrials+1:end) = pp(boolean(1-data_allT_perComb{iLoc, 1}));
            
            switch ip
                case 1, patch_both_perComb_p{iLoc} = pp_reorder;
                case 2, noise_both_perComb_p{iLoc} = pp_reorder;
                case 3, noiseMasked_both_perComb_p{iLoc} = pp_reorder;
            end
        end % end of iLoc
    end % end of ip
end % end of flagGetPatch

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% after this point, *_perBlk are no longer in use
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% [Behav] reorganize data so that PRS-trials are before the ABS-trials
data_both_allT_perComb = cell(nLoc8, nm7); % (1) itgtPrs, (2) answer, (3) correctness, (4) RT, (5) tgt phase

for iLoc = 1:nLoc8
    ii = data_allT_perComb{iLoc, 1}; % 1st layer is PrsInd (1 = tgt prs, 0 = tgt abs)
    
    % reorder data
    for im = 1 : nm7
        data_both_allT_perComb{iLoc, im} = [data_allT_perComb{iLoc, im}(boolean(ii)); data_allT_perComb{iLoc, im}(boolean(1-ii))];
    end
end

%% [Behav] derive SDT and pA based on ALL TRIALS
dprime_allT_perComb = nan(nLoc8, 1); % allT: calculated based on all trials
criterion_allT_perComb = dprime_allT_perComb;
pA_allT_perComb = nan(nLoc8, ntypes);
pC_allT_perComb = nan(nLoc8, ntypes);

for iLoc = 1:nLoc8
    if iLoc>5, switch iLoc, case 6, iLoc_all = [2,4]; case 7, iLoc_all = [3,5]; case 8, iLoc_all = 2:5; end
    else, iLoc_all = iLoc;
    end
    
    itgt_both = data_both_allT_perComb{iLoc, 1};  % prsInd
    resp_both = data_both_allT_perComb{iLoc, 2}; % answer (1=YES, 0=NO)
    correctness_both = data_both_allT_perComb{iLoc, 3}; % correctness (1=correct, 0=wrong)
    
    % pC
    pHit = mean(resp_both & itgt_both)*2; assert(pHit - mean(correctness_both & itgt_both)*2 < .000001)
    pFA = mean(resp_both & 1-itgt_both)*2; assert(pFA - mean((1-correctness_both) & (1-itgt_both))*2 < .000001)
    pC = mean(resp_both == itgt_both); assert(pC - mean(correctness_both) < .000001)
    pC_allT_perComb(iLoc, :) = [pHit, pFA, pC];
    
    % dprime, criterion
    [d,c] = SX_sim06_SDT(pHit, pFA, 1e5);
    dprime_allT_perComb(iLoc) = d;
    criterion_allT_perComb(iLoc) = c;
    
    % pA
    for itype = 1:ntypes
        respC_perComb = [];
        for ii = 1:length(iLoc_all)
            respC_perComb = [respC_perComb, squeeze(respC_perBlk(iCuedLoc_allT == iLoc_all(ii), itype, :))];
        end
        pA_allT_perComb(iLoc, itype) = nanmean(respC_perComb(:));
    end
end

%% [Behav] get ave/sem_perComb from _perBlk
cst_log_perBlk = log10(cst_perBlk);
RT_log_perBlk = log(RT_aveT_perBlk);

nm10 = 10;
data_perBlk_perComb = cell(nm10, nLoc8);
data_aveBlk_perComb = nan(nm10, nLoc8);
data_semBlk_perComb = data_aveBlk_perComb ;

for iLoc = 1:nLoc8
    if iLoc>5, switch iLoc, case 6, iLoc_all = [2,4]; case 7, iLoc_all = [3,5]; case 8, iLoc_all = 2:5; end
    else, iLoc_all = iLoc;
    end
    
    for im = 1:nm10
        switch im
            case 1, mm = dprime_perBlk;
            case 2, mm = criterion_perBlk;
            case 3, mm = RT_log_perBlk;
            case 4, mm = cst_log_perBlk;
            case 5, mm = pA_perBlk(3, :)'; % pA of all trials
            case 6, mm = pA_perBlk(1, :)'; % pA of PRS trials
            case 7, mm = pA_perBlk(2, :)'; % pA of ABS trials
            case 8, mm = pC_perBlk(3, :)'; % pC of ABS trials
            case 9, mm = pC_perBlk(1, :)';  % pHit
            case 10, mm = pC_perBlk(2, :)'; % pFA
        end
        dd = [];
        for ii = 1:length(iLoc_all)
            dd = [dd; mm(iCuedLoc_allT == iLoc_all(ii))];
        end
        data_perBlk_perComb{im, iLoc} = dd;
        data_aveBlk_perComb(im, iLoc) = mean(dd);
        data_semBlk_perComb(im, iLoc) = std(dd)/sqrt(nblocks);
    end
end

%% [Behav] RT median
% 1. PRS vs. ABS
% 2. correct vs. wrong
% 3. say yes vs. no

RT_log_allT_perLoc = cell(3, 2, nLoc5);
RT_log_mean_perLoc = nan(3, 2, nLoc5);

for im = 1:3
    for iLoc = 1:nLoc5
        RT_perLoc = squeeze(data_both_allT_perComb{iLoc, 4});
        ind1 = boolean(data_both_allT_perComb{iLoc, im} == 1); % PRS/correct/yes
        ind2 = boolean(data_both_allT_perComb{iLoc, im} == 0);% ABS/wrong/no
        % save RT of all trials
        RT_log_allT_perLoc{im, 1, iLoc} = RT_perLoc(ind1);
        RT_log_allT_perLoc{im, 2, iLoc} = RT_perLoc(ind2);
        % save the mean
        RT_log_mean_perLoc(im, 1, :) = mean(RT_perLoc(ind1));
        RT_log_mean_perLoc(im, 2, :) = mean(RT_perLoc(ind2));
    end
end

%% [patch-output] get CI
if flagGetPatch
    CI_mean = cell(1,nLoc8);
    CI_var = CI_mean;
    
    parfor iLoc = 1:nLoc8
        if iLoc>5, switch iLoc, case 6, iLoc_all = [2,4]; case 7, iLoc_all = [3,5]; case 8, iLoc_all = 2:5; end
        else, iLoc_all = iLoc;
        end
        noisePrs = noise_both_perComb_p{iLoc}(1:ntrials*length(iLoc_all));
        noiseAbs = noise_both_perComb_p{iLoc}(ntrials*length(iLoc_all)+1:end);
        
        idxHit = boolean(data_both_allT_perComb{iLoc, 2}(1:ntrials*length(iLoc_all))); % column 2 is the resp (1=yes, 0=no)
        idxFA = boolean(data_both_allT_perComb{iLoc, 2}(ntrials*length(iLoc_all)+1:end));
        idxMiss = boolean(1-idxHit);
        idxCR = boolean(1-idxFA);
        
        % hit (prs and say yes)
        noiseHit = cell2Mat3(noisePrs(idxHit));
        
        % FA (abs and say yes)
        if sum(idxFA) == ntrials, noiseFA = [];
        else, noiseFA = cell2Mat3(noiseAbs(idxFA));
        end
        
        % miss (prs and say no)
        if sum(idxMiss) == ntrials, noiseMiss = [];
        else, noiseMiss = cell2Mat3(noisePrs(idxMiss));
        end
        
        % CR (absent and say no)
        noiseCR = cell2Mat3(noiseAbs(idxCR));
        
        fprintf('%s: pHit=%.2f, pFA=%.2f, pMiss=%.2f, pCR=%.2f\n', namesLocComb{iLoc}, mean(idxHit), mean(idxFA), mean(idxMiss), mean(idxCR))
        
        % get CI
        CI_mean_ = mean(noiseHit,3) + mean(noiseFA,3) - mean(noiseCR,3) - mean(noiseMiss,3);
        CI_var_ = var(noiseHit,[],3) + var(noiseFA,[],3) - var(noiseCR,[],3) - var(noiseMiss,[],3);
        
        CI_mean{iLoc} = CI_mean_;
        CI_var{iLoc} = CI_var_;
    end % produce the 2 items defined above,
end

%% [patch-output] check power spectrum of patches
if flagGetPatch
    RC_checkNoiseFFT % output is fft1D_mean_both, 2x96 vector, should take 2:13 for SF = 1-4 cpd
end

%% save patch & related outputs
if flagGetPatch
    fprintf(' Saving patches...')
    save(nameTgtPatch, 'patch_both_perComb_p')
    save(namePatchOutputs, 'subjName', 'CI*', 'fft1D_mean_both') % search for x_record
    fprintf('  DONE\n')
end

%%
clear noise* patch* pp* nn*
clear subj filters* kernels* nB* nparams*  pV* R2*

%% save data
save(nameBehavMeas, '*_allT', 'n*', '*_perLoc', '*_perBlk', '*_perComb')
fprintf('  Data extracted and saved. \n')



