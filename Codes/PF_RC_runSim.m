
// clc,close all, clear all

// warning off
// format compact

// analysisMode = 1; % 1=simulation, 2=formal analysis
// addpath(genpath('fxn_analysis_RC'))
// addpath(genpath('fxn_exp'))
// addpath(genpath('fxn_simulation'))
// addpath(genpath('Data'))

// limit0to1 = @(x) min(max(x,0),1);

// %% set simulation params
// nsubj_all = 1;
// ntrials_all = [500, 1000]; % number of tgt-prs trials
// gaborCST_all = [.05, .1, .15, .2];

// ntrials_all = 1e3;%[1000]; % number of Gabor-PRS trials
// gaborCST_all = [1];

// nGaborCST = length(gaborCST_all);
// n_nsubj = length(nsubj_all);
// n_ntrials = length(ntrials_all);

// simPairs = combvec(nsubj_all,ntrials_all, gaborCST_all);
// nSimPairs = size(simPairs,2);

// %% step 1: get params
// SX_sim01_getParams
// noise.noiseCST = 0.2;
// % stim.gaborCST = 1;
// % stim.phase = 0; % phase of the target gabor; mute this line if want phase to be random on each trial
// noise.noiseProp = 0; % from AE2006, 0.035, 0.1, 0.19
// criterion_true = 0; % this is true after moving the two distributions by substracting their means
// iLoc = 0; % do not delete
// fprintf('noise cst = %.2f, noise prop. = %.2f, criterion = %d.\n', noise.noiseCST, noise.noiseProp, criterion_true)

// filtersSF_all = noise.filtersSF_all;
// noise.filtersSF_all = filtersSF_all;
// nfilters = length(filtersSF_all);
// filtersSF_all_log = log2(noise.filtersSF_all);
// fOri = 0:10:80;
// filtersOri_all = [-flip(fOri), 0, fOri] + 90;
// nfiltersOri = length(filtersOri_all);
// nfiltersSF = length(filtersSF_all);

// %% step 2: create filters at each SF channel
// plotFlag = 0;
// scr = [];
// [filter_sin, filter_cos] = SX_sim02_setFilters(stim, noise.filtersSF_all, filtersOri_all, plotFlag);
// disp('filters ready')

// %% make empty containers
// % simuContainer = cell(nsubj_all, n_ntrials, nGaborCST);
// % stims_both_allsubj = simuContainer;
// % stimPhase_allSubj = simuContainer;
// % energy_allSubj = simuContainer;
// % energy_norm_allSubj = simuContainer;
// % tgtPhase_allsubj = simuContainer;
// % behav_allSubj = simuContainer;
// % % noisePatch_allsubj = simuContainer;
// % CImean_allSubj = simuContainer;
// % CIvar_allSubj = simuContainer;
// % dprime_allSubj = simuContainer;
// % criterion_allSubj = simuContainer;
// % kernel_allSubj = simuContainer;
// % kernel_pred_allSubj = simuContainer;
// % params_est_allSubj = simuContainer;
// % nLL_allSubj = simuContainer;

// for ipair = 1:nSimPairs

//     nsubj = simPairs(1,ipair);
//     i_nsubj = find(nsubj == nsubj_all);

//     ntrials = simPairs(2,ipair);
//     i_ntrials = find(ntrials == ntrials_all);
//     ntrialsAll = ntrials*2;

//     stim.gaborCST = simPairs(3, ipair);
//     iTgtCST = find(stim.gaborCST == gaborCST_all);

//     fprintf('nsubj = %d (%d/%d), ntrials = %d (%d/%d), tgt cst = %.2f (%d/%d).\n', nsubj, i_nsubj, n_nsubj, ntrials,i_ntrials, n_ntrials,  stim.gaborCST, iTgtCST, nGaborCST)

//     for isubj = 1:nsubj
//         fprintf('  S%d/%d.', isubj, nsubj)
//         subjName = ['S',num2str(isubj)];

//         %% step 3: create stim
//         stim.phase = 0; %
//         [stims_both, tgtPhase, noisePatch_both] = SX_sim03_setStim(ntrials, stim, noise);
//         itarget = [ones(ntrials, 1); zeros(ntrials, 1)]; % present/absent at each trial

//         tgtPhase = [tgtPhase, zeros(1,ntrials)];

//         fprintf(' stimuli ready.')

//         stim_mean = nan(ntrialsAll, 1);
//         for in = 1:ntrialsAll
//             stim_mean(in) = mean(stims_both{in}(:));
//         end

//         %%
//         figure, hold on
//         histogram(stim_mean(1:ntrials), 'DisplayStyle', 'stairs','edgecolor', 'r', 'normalization', 'probability')
//         xline(mean(stim_mean(1:ntrials)), 'r');
//         histogram(stim_mean(ntrials+1:end), 'DisplayStyle', 'stairs', 'edgecolor', 'b', 'normalization', 'probability')
//         xline(mean(stim_mean(ntrials+1:end)), 'b');
//         %%
//         % figure, subplot(1,2,1), imshow(stims_both{1}), subplot(1,2,2), imshow(stims_both{ntrials+1})

//         %% step 4: compute energy
//         DV = nan(ntrialsAll,1);
//         energy = nan(ntrialsAll,nfiltersOri, nfiltersSF);
//         stimPhase = energy;

//         for itrial = 1:ntrialsAll
//             patch = stims_both{itrial};
//             % create the template
//             % Option 1: fixed phase at 0
//             stim.phase = 0;
//             % Option 2: random phase
//             %             stim.phase =  rand * 2 *pi;
//             % Option 3: phase from data
//             %             stim.phase =  tgtPhase(ss);
//             [template, ~] = exp_CreateGabor(stim, 1, 0);
//             template = limit0to1(template);

//             % get energy for RC analysis
//             DV(itrial) = patch(:)' * template(:);

//             for ifilterOri = 1:nfiltersOri
//                 for ifilterSF = 1:nfiltersSF
//                     [a, b] = SX_sim04_computeEnergy(patch, patch(1,1), filter_sin{ifilterOri, ifilterSF}, filter_cos{ifilterOri, ifilterSF}, 0, itrial <= ntrials, filtersSF_all(ifilterSF), filtersOri_all(ifilterOri));
//                     energy(itrial,ifilterOri, ifilterSF) = a;
//                     stimPhase(itrial,ifilterOri, ifilterSF) = b;
//                 end
//             end
//         end

//         energy_norm =  normEnergy(energy, ones(ntrialsAll,1) * stim.gaborCST);

//         %         energy_allSubj{isubj, i_ntrials, iTgtCST} = energy;
//         %         energy_norm_allSubj{isubj, i_ntrials, iTgtCST} = energy_norm;
//         %
//         fprintf(' energy computed.')


//         %% step 5: get behavior (for simulation only)
//         % way 0: Wyart etal., 2012
//         delta = .03; % additive bias
//         alpha = .1; % soft threshold
//         noiseInternal = 1.5; % internal noise SD
//         criterion = 16;

//         DV = (DV)/std(DV);
//         DV = DV+delta;
//         cstRespFxn = @(alpha) DV + alpha .* exp(-DV/alpha);
//         DV_nonL = cstRespFxn(alpha);
//         DV_nonL_noise = DV_nonL+randn(ntrialsAll, 1) * noiseInternal;

//         figure('Position', [0 300 600 200])
//         subplot(1,3,1), hold on
//         histogram(DV(1:ntrials)), histogram(DV(ntrials+1:end))
//         title('Energy norm + delta')

//         subplot(1,3,2), hold on
//         histogram(DV_nonL(1:ntrials)), histogram(DV_nonL(ntrials+1:end))
//         title('soft threshed')

//         subplot(1,3,3), hold on
//         histogram(DV_nonL_noise(1:ntrials)), histogram(DV_nonL_noise(ntrials+1:end))
//         xline(criterion, 'r-');
//         title('noise added')

//         resp = DV_nonL_noise>criterion;
//         mean(resp)

//         % way 1: eq 6 in AE2006
//         %         [template, ~] = exp_CreateGabor(stim, stim.gaborCST, 0);
//         %         DV_all = nan(1, ntrialsAll); % decision var
//         %         for ss = 1:ntrialsAll, DV_all(ss) = sum(stims_both{ss}(:) .* template(:)) + randn * noise.noiseProp; end
//         %         DV_all_norm = (DV_all - mean(DV_all))/std(DV_all);
//         %         figure, hold on, histogram(DV_all_norm(1:ntrials), 50, 'FaceColor', 'r', 'normalization', 'probability'), histogram(DV_all_norm(ntrials+1:end), 50, 'FaceColor', 'b', 'normalization', 'probability')
//         %         behav = DV_all_norm > criterion_true;

//         % way 2: RC
//         %         lumiBG = template(1,1);
//         %         energy_temp = nan(1,nfilters);
//         %         for ff = 1:nfilters
//         %             [a,b] = SX_sim04_computeEnergy(template , lumiBG, SFfilter_sin{ff}, SFfilter_cos{ff}, 0, ss<=ntrials, filterSF_all(ff));
//         %             energy_temp(ff) = a;
//         %         end
//         %         energy_temp_norm = energy_temp - mean(energy_temp);
//         %         SX_sim05_getBehavior

//         % way 3: original one
//         %         plotFlag = 0;
//         %         SX_sim05_getBehavior
//         %         behav = resp_noisy > criterion_true;

//         %         behav_allSubj{isubj, i_ntrials, iTgtCST} = resp;

//         % get CI
//         hitInd = resp & itarget;
//         FAInd = resp & ~itarget;
//         CRInd = ~resp & ~itarget;
//         missInd = ~resp & itarget;

//         %         ind4 = {hitInd; FAInd; CRInd; missInd};
//         %         NF_mean = cell(1,4); % NF: noise field
//         %         NF_var = cell(1,4);
//         %         for n = 1:4
//         %             pp = noisePatch_both(ind4{n});
//         %             NF_mean{n} = squeeze(nanmean(cat(3, pp{:}),3));
//         %             NF_var{n} = squeeze(nanstd(cat(3, pp{:}),[], 3));
//         %             if isempty(NF_mean{n}), NF_mean{n}=0;end
//         %             if isempty(NF_var{n}), NF_var{n}=0;end
//         %         end
//         %
//         %         CImean_allSubj{isubj, i_ntrials, iTgtCST} = NF_mean{1} + NF_mean{2} - (NF_mean{3} + NF_mean{4});
//         %         CIvar_allSubj{isubj, i_ntrials, iTgtCST} = NF_var{1} + NF_var{2} - (NF_var{3} + NF_var{4});

//         % step 6: SDT analysis
//         nHit = sum(hitInd);
//         nFA = sum(FAInd);
//         nCR = sum(CRInd);
//         pC = (nHit+nCR)/ntrialsAll;
//         [dprime, criterion] = SX_sim06_SDT(nHit, nFA, ntrials)

//         %         dprime2 = norminv(pC) * sqrt(2);

//         %         dprime_allSubj{isubj, i_ntrials, iTgtCST} = dprime;
//         %         criterion_allSubj{isubj, i_ntrials, iTgtCST} = criterion;

//         disp('Behavior simulated')

//         %% step 7: RC analysis
//         nB=1;

//         kernels2D_allB = nan(nB, ntypes, nfiltersOri, nfiltersSF); % 2D data
//         intercept_allB = kernels2D_allB;
//         R2_allB = nan(nB, ntypes, nfiltersOri, nfiltersSF);
//         R2_Tjur_allB = R2_allB;
//         pValues_allB = nan(nB, ntypes, nfiltersOri, nfiltersSF, 2); % 2 = p of intercept+p of slope
//         seprblity_kernel_allB = nan(nB, ntypes); % p values of check separability

//         for itype = 1:3   % PRS, ABS, BOTH
//             for iB = 1:nB
//                 % ========================================
//                 % decide the trial index
//                 if nB == 1 % NOT resample sample trials
//                     switch itype
//                         case 1, trialInd = 1:ntrials;
//                         case 2, trialInd =  (ntrials+1) : (ntrials*2);
//                         case 3, trialInd = 1:(ntrials*2);
//                     end
//                 else % resample sample trials
//                     switch itype
//                         case 1, trialInd = randi(ntrials, 1, ntrials);
//                         case 2, trialInd = randi(ntrials, 1, ntrials) + ntrials;
//                         case 3, trialInd = [randi(ntrials, 1, ntrials), randi(ntrials, 1, ntrials)+ntrials];
//                     end
//                 end
//                 % ========================================
//                 [kernel2D, ~, R2, R2_Tjur, pValues, intercept] = SX_sim07_RC(filtersSF_all, filtersOri_all, energy_norm(trialInd, :, :), resp(trialInd));

//                 % 1. save data
//                 intercept_allB(iB, itype, :, :) = intercept;
//                 kernels2D_allB(iB, itype, :, :) = kernel2D;
//                 R2_allB(iB, itype, :, :) = R2;
//                 R2_Tjur_allB(iB, itype, :, :) = R2_Tjur;
//                 pValues_allB(iB, itype, :, :, :) = pValues;

//                 % 2. separability
//                 seprblity_kernel_allB(iB, itype) = getSeparability(kernel2D);
//             end
//         end
//         kernels2D = squeeze(median(kernels2D_allB,1));
//         intercept = squeeze(median(intercept_allB,1));
//         R2 = squeeze(median(R2_allB,1));
//         R2_Tjur = squeeze(median(R2_Tjur_allB,1));
//         pValues = squeeze(median(pValues_allB,1));
//         seprblity_kernel = squeeze(median(seprblity_kernel_allB,1));

//         %
//         figure
//         for itype = 1:3
//             subplot(1,3,itype)
//             %             data2D_toPlot = flip(squeeze(data2D(iLoc , :, :))');
//             imagesc(filtersOri_all-90, flip(filtersSF_all_log), flip(squeeze(kernels2D(itype, :, :))'))

//             axis square
//             xline(0, 'r-', 'linewidth', 2);
//             yline(1, 'r-', 'linewidth', 2);
//             cc = colorbar;
//             %             caxis([caxisMin, caxisMax])

//             %             xticks(ticks_ORI)
//             %             xticklabels(ticklabels_ORI)
//             %             xlim(ticks_ORI([1,end]))
//             %
//             %             yticks(ticks_SF)
//             %             yticklabels(ticklabels_SF)
//             %             ylim(ticks_SF([1,end]))

//         end

//         %% marginalize
//         % fitMode = 2; % SSE = 1, MLE = 2;
//         % [kernel_pred, params_est, nLL] = SX_sim08_fit(filtersSF_all, kernel, model, fitMode);
//         % fprintf('estimated params: %.2f, %.2f, %.2f, %.2f\ndprime = %.2f, criterion = %.2f\n', params_est, dprime, criterion)

//         % kernel_allSubj{isubj, i_ntrials, iTgtCST} = kernel;
//         % kernel_pred_allPairs{i_nsubj, i_ntrials, iTgtCST} = kernel_pred_allSubj;
//         % params_est_allPairs{i_nsubj, i_ntrials, iTgtCST} = params_est_allSubj;
//         % nLL_allPairs{i_nsubj, i_ntrials, iTgtCST} = nLL_allSubj

//         % fprintf(' DONE!\n')
//     end
// end

// %% save data matrix
// % save(sprintf('fxn_simulation/dataMat_nCST%d_useNoise', noise.noiseCST*100))

// %% plot
// plot_sim

// makeBeep(ones(1,3) * .2, [1e3, 750, 500])