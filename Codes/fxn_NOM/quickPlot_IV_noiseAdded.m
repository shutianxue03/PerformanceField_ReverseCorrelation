
figure, hold on

h1 = histfit(IV_PRS); h1(1).FaceColor = 'r'; h1(1).LineWidth = .5; h1(1).FaceAlpha = .5; h1(2).Color = 'r'; h1(2).LineWidth = 1.5;

h2 = histfit(IV_ABS); h2(1).FaceColor = 'b'; h2(1).LineWidth = .5; h2(1).FaceAlpha = .5; h2(2).Color = 'b'; h2(2).LineWidth = 1.5;

h1 = histfit(IV_PRS_noiseAdd); 
h1(2).Color = 'k'; h1(2).LineWidth = 2; h1(2).LineStyle = '--'; 
h1(1).FaceColor = 'none'; h1(1).EdgeColor = 'none'; h1(1).LineWidth = 1; h1(1).FaceAlpha = .5; 

h2 = histfit(IV_ABS_noiseAdd); 
h2(2).Color = 'k'; h2(2).LineWidth = 2; h2(2).LineStyle = '--'; 
h2(1).FaceColor = 'none'; h2(1).EdgeColor = 'none'; h2(1).LineWidth = 1; h2(1).FaceAlpha = .5; 

xline(thresh, 'k-', 'linewidth', 2);
legend({'[PRS] IV', '[PRS] Gaussian fit', '[ABS] IV', '[ABS] Gaussian fit',...
    '[PRS] IV-noise added', '[ABS] Gaussian fit', '[ABS] IV-noise added', '[ABS] Gaussian fit', ...
    'Threshold'})