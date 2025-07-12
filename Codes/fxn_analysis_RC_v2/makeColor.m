
figure
getFlag = 0;
while  ~getFlag
    colors = rand(3,nsubj);
    hold on
    for isubj = 1:nsubj
        plot(isubj, 1, 'o', 'MarkerFaceColor', colors(:, isubj), 'MarkerEdgeColor', 'w', 'MarkerSize', 10)
    end
    xlim([.5, nsubj+.5])
    getFlag = input('accept this color? (1=YES, 0=NO): ');
end

save(fileName_color, 'colors')

close all
clc