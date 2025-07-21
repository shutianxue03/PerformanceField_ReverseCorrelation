
% close all
nLines = 2; % there are two lines in each location combination
figure('Position', [0 0 2000 800])

for ifeature = 1:nFeatures
    switch ifeature
        case 1, marg = margORI; margPred = margPred_ORI; params = margParams_ORI;
        case 2, marg = margSF; margPred = margPred_SF; params = margParams_SF;
    end

    for iType = 1:nTypes

        subplot(nFeatures, nTypes, (ifeature-1)*nTypes+iType), hold on

        if flag_standEnergy
            if ifeature==1, ylim([-.05, .2]), xlim([-90, 90]), xticks(-90:45:90), title(namesType{iType}), xline(0, 'color', ones(1,3)*.5);
            else, ylim([-.05, .1]), xlim([0,2]), xline(1, 'color', ones(1,3)*.5); xticks(0:.5:2), xticklabels(round(2.^(0:.5:2), 2)),
            end
        else
            if ifeature==1, ylim([-3, 7]), xlim([-90, 90]), xticks(-90:45:90), title(namesType{iType}), xline(0, 'color', ones(1,3)*.5);
            else, ylim([-1, 2]), xlim([0,2]), xline(1, 'color', ones(1,3)*.5); xticks(0:.5:2), xticklabels(round(2.^(0:.5:2), 2)),
            end
        end

        yline(0, 'color', ones(1,3)*.5);

        for iLine = 1:nLines
            color = colors_comb(iLocComb_all(iLine), :);
            if flag_interpolate>1
                axis_ln = axis_tuning{ifeature};
                if ifeature==1 % only interpolate the center
                    nNotItp = sum(axis_ln<=-30);
                    axis_ln_itp = linspace(axis_ln(nNotItp), axis_ln(end-nNotItp+1), flag_interpolate);
                    x = [axis_ln(1:nNotItp-1), axis_ln_itp, axis_ln(end-nNotItp+2:end)];
                else
                    x = linspace(axis_ln(1), axis_ln(end), flag_interpolate);
                end
            else
                x=axis_tuning{ifeature};
            end

            % plot the data
            plot(x, squeeze(marg(iLine, iType, :)), 'o', 'color', color)
            % plot the pred
            plot(x, squeeze(margPred(iLine, iType, :)), '-', 'color', color)

            % plot SF peak
            if ifeature==2
                switch ifamily_perF(2)
                    case 2
                        xline(log2(margParams_SF(iLine, iType, 1)), 'color', color);
                    case 12
                        xline(log2(margParams_SF(iLine, iType, 1)), '-', 'color', color);
                        xline(log2(margParams_SF(iLine, iType, 5)), '--', 'color', color);
                end
            end
        end % iLine
    end % iType
    set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)
end % ifeature
sgtitle(sprintf('%s iB=%d', subjName, iB))

