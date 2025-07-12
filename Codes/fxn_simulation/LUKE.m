%% 0: INITIALIZE %%
close all; clear all; clc;
cols = cbrewer('qual', 'Set1', 4);

%% 1: SETUP STIMULI 

% Stimulus variables
sf_range = logspace(log10(1),log10(3.5),39); % Range of SFs used to compute energy
ori_range = -90:5:90;                        % Range of orientations used
n_ori = length(ori_range);                   % Number of orientation bins
n_sf = length(sf_range);                      % Number of SF bins
% sf_mdpt = find(round(sf_range,2) == 2);      % Index where SF == 2
ori_mdpt = ceil(n_ori/2);                   % Index where orientation == 0

% Using stimuli generated from my experiment code.
% 20% contrast vertical 2cpd Gabors, 3 degrees of visual angle in size, embedded
% in white noise band-pass filtered to contain SFs between 0.5 and 4 cpd.
% Orientation content was not filtered.
% 3000 present trials, 3000 absent trials.
load('presentStims.mat');
load('absentStims.mat');

% Plot example stimuli
figure('Position', [0 0 1200 600]); hold on;
subplot(1,2,1); hold on;
imagesc(1:96,1:96,squeeze(presentStims(5,:,:)))
xlim([1 96])
ylim([1 96])
colormap('gray')
caxis([0.25 0.8])
title('Example present stimulus')
set(gca,'XColor', 'none','YColor','none')
set(gca, 'FontSize', 20)
subplot(1,2,2); hold on;
imagesc(1:96,1:96,squeeze(absentStims(5,:,:)))
xlim([1 96])
ylim([1 96])
colormap('gray')
title('Example absent stimulus')
caxis([0.25 0.8])
set(gca,'XColor', 'none','YColor','none')
set(gca, 'FontSize', 20)


% Only want to do reverse-correlation on orientation so things are simpler
% and I don't have to deal with 2D plotting and programming. 
% I take the mean across SF energy, rather than just a slice, because it
% makes things less noisy, and helps with the fact that my noise is
% filtered to have the highest power below 2cpd.
energy_prs = squeeze(mean(energy_prs,2));
energy_abs = squeeze(mean(energy_abs,2));
ntrial = size(energy_prs,2)*2;               % Number of total trials
energy = cat(2,energy_prs,energy_abs);  % Put stimuli together

% For the beta-weight estimation, we first want to subtract the signal
% (Gabor stimulus) out of the energy. This way we isolate the effects of
% fluctuations in noise on behavior, independent of the signal. We take the
% mean of each cell in the energy matrix (i.e., each orientation),
% and subtract the mean across trials.
% Second, we set the standard deviation of each cell to 1 across trials so
% that our beta-weight estimates are comparable across the kernel.
energy_prs_norm = (energy_prs - mean(energy_prs,2))./ std(energy_prs,[],2);
energy_abs_norm = (energy_abs - mean(energy_abs,2))./ std(energy_abs,[],2);
energy_norm = cat(2,energy_prs_norm,energy_abs_norm);              % Put stimuli together

%% 2: SETUP OBSERVER
% Make an 'attention filter' function. 
% It creates a 1D gaussian PDF across orientations for -90 to 90 degrees,
% which we will use to compute a response to the stimulus energy.
% We're assuming that attention only changes the bandwidth, so we'll allow
% for changing the standard deviation of the gaussian across each
% dimension, but assume that the mean is centered on the stimulus
% orientation (0 degrees)
makeAttnFilt = @(sigmaORI) normpdf(ori_range, ori_range(ori_mdpt), sigmaORI);

% I normalize each distribution so the max is equal to 1. This way the only
% thing that differs between cue conditions is the width, not the magnitude
% of response at the peak.

% Filter for neutral trials
baselineSigmaORI = 30;                        % Distributed attention ORI tuning as a baseline
neutralFilt = makeAttnFilt(baselineSigmaORI); % Get gaussian filter
neutralFilt = neutralFilt./max(neutralFilt);  % Normalize peak to 1

% Filter for valid trials
% Half as much spread as neutral:
validFilt = makeAttnFilt(baselineSigmaORI*0.5);
validFilt = validFilt./max(validFilt); 

% Filter for invalid trials
% Twice as much spread for invalid:
invalidFilt = makeAttnFilt(baselineSigmaORI*2);
invalidFilt = invalidFilt./max(invalidFilt); 

% Internal noise level.
% I got this value by hand - it only changes the magnitude of the
% beta-weights and d', but doesn't change trends.
intNoiseLvl = 0.08;


%% 3: RUN TRIALS %%
% Point-wise multiply each stimulus with the attention filter, then sum
% across all orientations to get a point estimate.
% Then, we add random noise with a mean of 0 to each point estimate.
resp_Neutral = sum(energy .* neutralFilt',1) + randn(1,ntrial)*intNoiseLvl;  % Neutral
resp_Valid = sum(energy .* validFilt',1) + randn(1,ntrial)*intNoiseLvl;      % Valid
resp_Invalid = sum(energy .* invalidFilt',1) + randn(1,ntrial)*intNoiseLvl;  % Invalid

% The different gaussian filters result in different magnitudes of response
% estimates. To properly equate things across cue conditions, I subtract
% the mean out of each response. Then, the criterion can be set to 0 for
% all conditions.
resp_Neutral = resp_Neutral - mean(resp_Neutral);
resp_Valid = resp_Valid - mean(resp_Valid);
resp_Invalid = resp_Invalid - mean(resp_Invalid);
criterion = 0;

% Get reponses; > 1 = respond present, <= 1 = response absent
behav_Neutral = resp_Neutral > criterion;
behav_Valid = resp_Valid > criterion;
behav_Invalid = resp_Invalid > criterion;

%% 4. SDT ANALYSES %%

% Get hit-rate and false-alarm rate for each cue condition
neutralSDT = [mean(behav_Neutral(1:(ntrial/2))) mean(behav_Neutral(((ntrial/2)+1):end))];
validSDT = [mean(behav_Valid(1:(ntrial/2))) mean(behav_Valid(((ntrial/2)+1):end))];
invalidSDT = [mean(behav_Invalid(1:(ntrial/2))) mean(behav_Invalid(((ntrial/2)+1):end))];

% Get d' for each condition using z-transform
dprime_Neutral = norminv(mean(behav_Neutral(1:(ntrial/2)))) - norminv(mean(behav_Neutral(((ntrial/2)+1):end)));
dprime_valid = norminv(mean(behav_Valid(1:(ntrial/2)))) - norminv(mean(behav_Valid(((ntrial/2)+1):end)));
dprime_invalid = norminv(mean(behav_Invalid(1:(ntrial/2)))) - norminv(mean(behav_Invalid(((ntrial/2)+1):end)));

% Plot d'
figure; hold on;
LH_scatter(1, dprime_valid, cols(2,:), 500)
LH_scatter(2, dprime_invalid, cols(1,:), 500)
LH_scatter(3, dprime_Neutral, 'k', 500)
xticks(1:3)
xticklabels({'Valid', 'Invalid', 'Neutral'})
xlim([0.5 3.5])
ylabel('d prime')
set(gca,'FontSize',30)

%% 5. RC ANALYSES %%
% Initialize collection arrays for beta-weight and error (SE on beta-weight
% estimate) for each cue condition.
kernel_neutral = nan(1,n_ori);
kernel_valid = nan(1,n_ori);
kernel_invalid = nan(1,n_ori);

err_neutral = nan(1,n_ori);
err_valid = nan(1,n_ori);
err_invalid = nan(1,n_ori);

% Loop through each orientation and compute beta-weight and error
% Here, I'm folding across orientations (i.e., treating -20 and 20 degrees
% as the same) because we know the simulated observer treats these the
% same, and it will make the data less noisey.
for i_ori = 1:n_ori
    
    % Valid
    [b,~,stat] = glmfit([energy_norm(i_ori,:)'; energy_norm((n_ori+1)-i_ori,:)'] ,[behav_Valid'; behav_Valid'],'binomial','link','logit');
    kernel_valid(i_ori) = b(2);
    err_valid(i_ori) = stat.se(2);
    
    % Invalid
    [b,~,stat] = glmfit([energy_norm(i_ori,:)'; energy_norm((n_ori+1)-i_ori,:)'] ,[behav_Invalid'; behav_Invalid'],'binomial','link','logit');
    kernel_invalid(i_ori) = b(2);
    err_invalid(i_ori) = stat.se(2);
    
    % Neutral
    [b,~,stat] = glmfit([energy_norm(i_ori,:)'; energy_norm((n_ori+1)-i_ori,:)'] , [behav_Neutral'; behav_Neutral'],'binomial','link','logit');
    kernel_neutral(i_ori) = b(2);
    err_neutral(i_ori) = stat.se(2);
    
end

%% 6. RC PLOT %%

% To make the plotting nicer, I spline interpolate the beta-weight and
% error estimates to have a much finer-grain (1000 points rather than 39)
oriRangeInterp = linspace(-90,90,1000);
validKernelInterp = interp1(ori_range,kernel_valid,oriRangeInterp,'Spline');
invalidKernelInterp = interp1(ori_range,kernel_invalid,oriRangeInterp,'Spline');
neutralKernelInterp = interp1(ori_range,kernel_neutral,oriRangeInterp,'Spline');
validErrInterp = interp1(ori_range,err_valid,oriRangeInterp,'Spline');
invalidErrInterp = interp1(ori_range,err_invalid,oriRangeInterp,'Spline');
neutralErrInterp = interp1(ori_range,err_neutral,oriRangeInterp,'Spline');

% Plot
figure; hold on;
plot([ori_range(1) ori_range(end)], [0 0], 'LineWidth', 2, 'color', [0.5 0.5 0.5], 'Handlevisibility', 'off')
plot([0 0], [-2 2], 'LineWidth', 2, 'color', [0.5 0.5 0.5], 'Handlevisibility', 'off')
plot(oriRangeInterp, validKernelInterp, 'color', cols(2,:),'LineWidth',3)
plot(oriRangeInterp, invalidKernelInterp, 'color', cols(1,:),'LineWidth',3)
plot(oriRangeInterp, neutralKernelInterp, 'color', 'k','LineWidth',3)
legend('Valid', 'Invalid', 'Neutral')
patch([oriRangeInterp fliplr(oriRangeInterp)], [validKernelInterp + validErrInterp,...
    fliplr(validKernelInterp - validErrInterp)],...
    cols(2,:),'EdgeColor', 'none', 'FaceAlpha', .2);
patch([oriRangeInterp fliplr(oriRangeInterp)], [invalidKernelInterp + invalidErrInterp,...
    fliplr(invalidKernelInterp - invalidErrInterp)],...
    cols(1,:),'EdgeColor', 'none', 'FaceAlpha', .2);
patch([oriRangeInterp fliplr(oriRangeInterp)], [neutralKernelInterp + neutralErrInterp,...
    fliplr(neutralKernelInterp - neutralErrInterp)],...
    'k','EdgeColor', 'none', 'FaceAlpha', .2);
LH_scatter(ori_range,kernel_valid,cols(2,:))
LH_scatter(ori_range,kernel_invalid,cols(1,:))
LH_scatter(ori_range,kernel_neutral,'k')
xlim([ori_range(1) ori_range(end)])
xlabel('Orientation')
ylabel('Sensitivity (a.u.)')
ylim([-0.1 1.3])
set(gca,'FontSize',20)
yticks([0 0.5 1])
xticks([-80:40:80])


