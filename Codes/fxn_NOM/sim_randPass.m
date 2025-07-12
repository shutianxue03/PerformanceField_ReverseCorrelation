% Repeat the random sampling of trials to create the double pass and generate the distribution of pA_PRS/ABS
%
close all

ntrials_sim = 1e4;

% mu_PRS = 1;
% sigma_PRS = 1;
% mu_ABS = -1;
% sigma_ABS = 1;
% thresh = 0;

mu_all = -1:1;
sigma_all = .5:.5:1.5;
thresh = -.5:.5:.5;

set_all = combvec(mu_all, sigma_all, thresh);
ns = size(set_all, 2);
for is = 1:ns
    mu_PRS = set_all(1, is);
    sigma_PRS = set_all(2, is);
    thresh = set_all(3, is);
    
    IV_PRS = randn(1,ntrials_sim) * sigma_PRS + mu_PRS;
    % IV_ABS = randn(1,ntrials_sim) * sigma_ABS + mu_ABS;
    
    IR_PRS = IV_PRS >= thresh; % so that 1 is the correct answer
    % IR_ABS = IV_ABS >= thresh; % so that 0 is the correct answer
    
    IR_PRS_pass = reshape(IR_PRS, [2, ntrials_sim/2]);
    % IR_ABS_pass = reshape(IR_ABS, [2, ntrials_sim/2]);
    
    pA_PRS = mean(IR_PRS_pass(1,:) == IR_PRS_pass(2,:));
    % pA_ABS = mean(IR_ABS_pass(1,:) == IR_ABS_pass(2,:));
    
    %%
    ni = 1e4;
    pA_PRS_all = nan(1, ni);
    % pA_ABS_all = pA_PRS_all;
    
    parfor ii=1:ni
        pA_PRS_ = mean(IR_PRS(randperm(ntrials_sim, ntrials_sim/2)) == IR_PRS(randperm(ntrials_sim, ntrials_sim/2)));
        %     pA_ABS_ = mean(IR_ABS(randperm(ntrials_sim, ntrials_sim/2)) == IR_ABS(randperm(ntrials_sim, ntrials_sim/2)));
        
        pA_PRS_all(ii) = pA_PRS_;
        %     pA_ABS_all(ii) = pA_ABS_;
    end
    
    figure('Position', [2000 500 1000 400])
    subplot(2,1,1)
    histogram(IV_PRS)
    xline(thresh, 'r-', 'linewidth', 2);
    xlim([-3, 3])
    title(sprintf('mu=%.2f sigma=%.2f thresh=%.2f', mu_PRS, sigma_PRS, thresh))

    subplot(2,1,2), hold on
    histogram(pA_PRS_all)
    xlim([.5, 1])
    xline(pA_PRS, 'r-', 'linewidth', 2);
    set(findall(gcf, '-property', 'fontsize'), 'fontsize',18)
    mu_PRS
    % title(sprintf('%.2f (%.4f)', median(pA_PRS_all), std(pA_PRS_all)))
    % subplot(1,2,2), histogram(pA_ABS_all), xline(pA_ABS, 'linewidth', 2); title(sprintf('%.2f (%.4f)', median(pA_ABS_all), std(pA_ABS_all)))
end