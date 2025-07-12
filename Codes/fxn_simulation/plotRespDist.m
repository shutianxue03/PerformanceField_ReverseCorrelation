figure('Position',[0 300 400 1200])
suptitle(locName,8)
ppd = 32;

%%
subplot(4,1,1),hold on
plot(stepSF,(energy(1:ntrials/2,:))') % present trials
plot(stepSF,(energy(ntrials/2+1:end,:))','k-') % absent trials are black lines
plot([target_SF,target_SF],ylim,'k-')
xlabel('SF channels (cpd)'),ylabel('Energy')
title('Un-normalized energy')

%%
subplot(4,1,2), hold on
plot(stepSF,tuning,'o-'),
plot([target_SF,target_SF],ylim,'k-')
xlabel('SF channels (cpd)'),ylabel('Sensitivity')
title('Internal tuning function')

%%
subplot(4,1,3),hold on
h_prs = histogram(resp(1:ntrials/2),ntrials/10,'Normalization','Probability');
h_abs = histogram(resp(ntrials/2+1:end),ntrials/20,'Normalization','Probability');
plot([0 0], ylim,'k-','LineWidth',2)
legend('Present trials','Absent trials','Criterion (response = 0)')
xlabel('Response'),ylabel('Proportion')
title('Distribution of prs/abs trials')

%%
subplot(4,1,4)
hold on
plot(stepSF,(energy_norm(1:ntrials/2,:))') % present trials
plot(stepSF,(energy_norm(ntrials/2+1:end,:))') % absent trials are black lines
plot([target_SF,target_SF],ylim,'k-')
xlabel('SF channels'),ylabel('Energy')
title('Normalized energy')
