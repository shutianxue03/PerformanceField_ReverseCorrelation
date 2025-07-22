
function make_sgtitle(title_, subjName, nB, nsubj, nAllTrials)

ntrials_mean = round(mean(nAllTrials));
ntrials_SEM = round(std(nAllTrials)/sqrt(nsubj));

if nsubj == 1
    sgtitle(sprintf('%s\n%s: %d trials per loc (nB=%d)', title_, subjName, nAllTrials, nB), 'FontSize',22)
else
    sgtitle(sprintf('%s\nn = %d (%d +- %d) (nB=%d)', title_, nsubj, ntrials_mean, ntrials_SEM, nB), 'FontSize',22)
end