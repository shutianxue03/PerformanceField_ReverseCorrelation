

flag_plotDiff = 1;

sz_fig = [400, 400];


if ifeature== 1
    switch iLocComb_all(1)
        case  6, y_ticks = linspace(0, 1.6, 5)/1e4;
        case  5, y_ticks = linspace(0, 4, 5)/1e4;
    end
    y_ticks = linspace(0, 4, 5)/1e4;
else
    switch iLocComb_all(1)
        case 6, y_ticks = linspace(0, 1.6, 5)/1e4;
        case 5, y_ticks = linspace(0, 3.6, 5)/1e4;
    end
    y_ticks = linspace(0, 3.6, 5)/1e4;
end

for iIC = 1:nIC_ % if using CV, nIC_=1
    dev_allF = nan(length(ifamily_all), nsubj);
    for ifamily = ifamily_all
        iif = find(ifamily == ifamily_all);
        dev_allF(iif, :) = dev_bestGroup_perModePerFamily{ifamily, MCmode, iIC};
    end
    
    % get delta
    dev_allF = dev_allF - min(dev_allF(:));
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    basicFxn_drawBars(dev_allF', [], [0 0 0;0,0,0], namesFamily_all(ifamily_all), y_ticks, [], flag_plotIDVD, flag_plotDiff, ...
        sprintf('[%s] %s L%d%d [ORI%d SF%d]', namesMCmode{MCmode}, namesFeature{ifeature}, iLocComb_all, nORI, nSF), 0, sz_fig);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% save figure
    switch length(ifamily_all)
        case 2, saveas(gcf, sprintf('%scompFamily_n%d_M%d%d_%s.jpg', nameFolder_fig, nsubj, ifamily_all, namesMCmode{MCmode}))
        case 3, saveas(gcf, sprintf('%scompFamily_n%d_M%d%d%d_%s.jpg', nameFolder_fig, nsubj, ifamily_all, namesMCmode{MCmode}))
    end
    
end % iIC





