%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot6_params_grant.m
% Last updated by Shutian Xue on 07/22/2025
%
% Description:
%   This script plots different tuning characteristics as a function of locations 
%   and computes group statistics and effect sizes for each bootstrap.
%
% Usage:
%   Called from plotAll_part1.m after compiling bootstrapped and fitted parameter data.
%
% Inputs (from workspace):
%   - margTuningC_ORI_allSubj, margTuningC_SF_allSubj: tuning characteristics for all subjects
%   - namesFeature, namesTunC_unit_perF, ifamily_perF, paramMode, nFeatures, nsubj, etc.
%
% Outputs:
%   - Saves summary bar plots of parameter distributions to the
% 
% Helper functions:
%   - basicFxn_drawBars (function)
%   - getCI (function)

% Plotting settings
colors = colors_comb(iLocComb_all, :);
x_ticks = namesLocComb(iLocComb_all);
paramMode = 2;
sz_fig = [350 350]; % size of the figure canvas

% Define the directory to save figures
nameFolder_Fig_tunC = sprintf('%s/%s/params/%s/', nameFigFolder, name_numFilters_Fitting, nameFileLoc);
if isempty(dir(nameFolder_Fig_tunC)), mkdir(nameFolder_Fig_tunC), end

% Loop through each feature
for iFeature = 1:nFeatures

    % Define the x-axis and parameters for the current feature
    xaxis = axis_tuning{iFeature};
    nfilters = length(xaxis);
    ifamily = ifamily_perF(iFeature);
    namesTunC = namesTunC_unit_perF{ifamily, paramMode};
    nTunC_full = length(namesTunC);

    y_ticks_all = [];

    %% Obtain y_ticks of ORI/SF domain (3 values)
    switch ifamily_perF(iFeature)
        case 1
            if flag_plotIDVD, y_ticks_all{1} = linspace(0, 36, 5); % preferred ori
            else, y_ticks_all{1} = linspace(0, 20, 5); % preferred ori
            end

            if flag_plotIDVD
                y_ticks_all{2} = linspace(0, .2, 5); % peak amp
                y_ticks_all{3} = linspace(12, 56, 5); % band
                y_ticks_all{4} = linspace(-.05, .03, 5); % baseline
            else
                y_ticks_all{2} = linspace(0, .12, 5); % peak amp
                y_ticks_all{3} = linspace(10, 50, 5); % band
                y_ticks_all{4} = linspace(-.05, .03, 5); % baseline
            end
        case 8
            % for L18
            y_ticks_all{1} = linspace(0, 16, 5); % preferred ori
            y_ticks_all{2} = linspace(.04, .12, 5); % peak amp
            y_ticks_all{3} =  linspace(30, 90, 5); % trough ori
            y_ticks_all{4} =  linspace(-.02, .02, 5); % trough amp
            y_ticks_all{5} = linspace(30, 50, 5); % band
            y_ticks_all{6} = linspace(-.02, .02, 5); % baseline

            % for others
            y_ticks_all{1} = linspace(0, 16, 5); % preferred ori
            y_ticks_all{2} = linspace(0, .12, 5); % peak amp
            y_ticks_all{3} =  linspace(20, 60, 5); % trough ori
            y_ticks_all{4} =  linspace(-.02, .02, 5); % trough amp
            y_ticks_all{5} = linspace(20, 60, 5); % band
            y_ticks_all{6} = linspace(-.02, .02, 5); % baseline

        case 2 % log parabola
            if flag_plotIDVD
                y_ticks_all{1} = linspace(0, 2, 5); % peak SF
                y_ticks_all{2} = linspace(.01, .09, 5); % peak amp
                y_ticks_all{3} = linspace(0, 2, 5); % bandwidth
                y_ticks_all{4} = linspace(-.08, 0, 5); % baseline
            else
                y_ticks_all{1} = linspace(0, 2, 5); % peak SF
                y_ticks_all{2} = linspace(.02, .06, 5); % peak amp
                y_ticks_all{3} = linspace(.3, 1.5, 5); % bandwidth
                y_ticks_all{4} = linspace(-.08, 0, 5); % baseline
            end
        case 3 % truncated log parabola
            y_ticks_all{1} = linspace(0, 2, 5); % peak SF
            y_ticks_all{2} = 0:.02:.08; % peak amp
            y_ticks_all{3} = .3:.2:1.1; % bandwidth
            y_ticks_all{4} = linspace(-.04, .04, 5); % baseline
            y_ticks_all{5} = 0:.03:.12; % trunc

        case 12 % double peak
            % if paramMode==ticks{2} = 0:.12:.24; else, ticks{2} = 0:.06:.12; end % gain vs. peak amp
            % if flag_plotOctave, ticks{3} = .4:1.2:2.8; else, ticks{3} = 1:.5:3.5; end % octave vs. cpd
            y_ticks_all{1} = linspace(0, 2, 5); % peak SF 1
            y_ticks_all{2} = 0:.08:.16; % peak amp 1
            y_ticks_all{3} = .3:.4:1.1; % bandwidth 1
            y_ticks_all{4} = [.678, 1.339, 2]; % peak SF2 (on log2 scale)
            y_ticks_all{5} = 0:.06:.12; % peak amp 2
            y_ticks_all{6} = .1:.5:1.1; % bandwidth 2

    end

    %% Loop through each tuning characteristic
    for iTunC = 1:nTunC_full

        text_title = sprintf('%s %s', namesFeature{iFeature}, namesTunC{iTunC});
        
        % Compute median for each observer
        switch iFeature
            case 1, data_allSubj = margTuningC_ORI_allSubj(:, :, :, itype, iTunC); % nsubj x nB x nLoc2 x ntypes x nTunC
            case 2, data_allSubj = margTuningC_SF_allSubj(:, :, :, itype, iTunC);
        end
        med_allSubj = getCI(data_allSubj, 1, 2);

        % Obtain ytick labels and ref
        y_ticklabels = y_ticks_all{iTunC};
        ref = nan;
        switch ifamily_perF(iFeature)
            case 1, if find(iTunC==[1,4]), ref = 0; end
            case 8
                switch iTunC
                    case 1, ref = 0;
                    case 6, ref = 0;
                end

            case 2 % log parabola
                switch iTunC
                    case 1, y_ticklabels = round(2.^y_ticks_all{iTunC}, 1); ref = log2(2); med_allSubj = log2(med_allSubj); % pref SF
                    case 4, ref = 0; % baseline
                end
            case 3 % truncated log parabola
                switch iTunC
                    case 1, y_ticklabels = round(2.^y_ticks_all{iTunC},2); ref = log2(2); med_allSubj = log2(med_allSubj); % pref SF
                    case 4, ref = 0; % baseline
                end
            case 12 % double peak
                if find(iTunC==[1,4]), y_ticklabels = round(2.^y_ticks_all{iTunC},2); ref = log2(2);  end

        end

        % Obtain CI of t-test statistics (to print in the title)
        for iB = 1:nB
            x=squeeze(data_allSubj(:, iB, :));[~, p, ~, stats] = ttest(x(:, 1), x(:, 2));cohenD = fxn_getES(x(:, 1), x(:, 2));
            p_allB(iB) = p; t_allB(iB) = stats.tstat; cohenD_allB(iB) = cohenD;
        end
        [p_med, p_lb, p_ub] = getCI(p_allB, 1, 2);
        [t_med, t_lb, t_ub] = getCI(t_allB, 1, 2);
        [d_med, d_lb, d_ub] = getCI(cohenD_allB, 1, 2);

        % Define the title for the plot
        title_CI = sprintf('t=%.3f [%.3f, %.3f], p=%.3f [%.3f, %.3f], d=%.3f [%.3f, %.3f]', t_med, t_lb, t_ub, p_med, p_lb, p_ub, d_med, d_lb, d_ub);

        % Plot the data
        text_title = sprintf('%s\n%s', text_title, title_CI);
        %------------------------------%
        flag_sig = basicFxn_drawBars(med_allSubj, ref, colors, x_ticks, y_ticks_all{iTunC}, y_ticklabels, flag_plotIDVD, flag_plotDiff, text_title, 0, sz_fig);
        %------------------------------%

        % Save the figure
        saveas(gcf, sprintf('%sn%d_%s%d%s.jpg', nameFolder_Fig_tunC, nsubj, namesFeature{iFeature}, iTunC, flag_sig))

    end % end of iparam
end