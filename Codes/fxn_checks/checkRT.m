
%% check if RT is an indicator of confidence
iLoc = 1;
RT_both_log = data_both(:, iLoc, 4);
correct_both = data_both(:, iLoc, 3);

figure, hold on
histogram(RT_both_log(1:ntrials))
histogram(RT_both_log(ntrials+1:end))
xticks(-4:0)
xticklabels(round(exp(-4:0)*1000))
xlim([-4,0])
xlabel('RT (ms)')

RT_log_edge = quantile(RT_both_log(1:ntrials), 0:.25:1);
pC = nan(length(RT_log_edge), 1);
for iedge = 1:length(RT_log_edge)-1
    pC(iedge) = mean(correct_both(RT_both_log <= RT_log_edge(iedge+1) & RT_both_log >= RT_log_edge(iedge)));
end

figure
bar(pC)
% makes sense, longer RT, lower accuracy
