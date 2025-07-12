% clc

if size(threshT_allSess, 1) < 2, error('ALERT: no qstair data!!!'), end
lastThreshUsed = threshT_allSess(end-1, :);
currentThreshMeasured = threshT_allSess(end, :);
% folderName = sprintf('%s/Data',pwd);
disp('==============')
disp('** last thresh used **')
fprintf('             %.4f\n%.4f %.4f %.4f \n             %.4f\n\n', lastThreshUsed([3,2,1,4,5]))
disp(lastThreshUsed)

disp('** current thresh measured **')
fprintf('             %.4f\n%.4f %.4f %.4f \n             %.4f\n\n', currentThreshMeasured([3,2,1,4,5]))
disp(currentThreshMeasured)


comp = (currentThreshMeasured-lastThreshUsed)./lastThreshUsed;
fprintf('** compare threshold **\n             %.4f\n%.4f %.4f %.4f \n             %.4f\n\n', comp([3,2,1,4,5]))


