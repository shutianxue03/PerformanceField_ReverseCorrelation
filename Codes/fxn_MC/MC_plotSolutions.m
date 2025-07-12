
% plot all the solutions
ns = length(solution);
params_allRep = nan(ns, nparamsAll); 
RSS_allRep = nan(ns,1);
for i = 1:ns
    params_allRep(i, :) = solution(i).X;
    RSS_allRep(i) = solution(i).Fval;
end

figure('Position', [0 200 nparamsAll*300 200])
for iparam = 1:nparamsAll
    subplot(1,nparamsAll+1, iparam), hold on
    bar(params_allRep(:, iparam), 'EdgeColor', 'w')
    title(sprintf('Parameter #%d', iparam))
    bar(params_allRep(1, iparam), 'FaceColor', 'r')
end

% hist of error
subplot(1,nparamsAll+1, nparamsAll+1), hold on
bar(RSS_allRep, 'EdgeColor', 'w')
bar(RSS_allRep(1), 'FaceColor', 'r')
title('RSS')
