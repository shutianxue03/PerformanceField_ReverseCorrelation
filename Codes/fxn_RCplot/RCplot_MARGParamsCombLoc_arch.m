
% close all
clear ydiffErr

% empty containers
params_all = cell(2, ntypes, max([nparamsOri, nparamsSF])); 
params_med_all = params_all;
params_SEM_neg_all = params_all;
params_SEM_pos_all = params_all;
params_groupAVE_all = params_all;
params_groupSEM_all = params_all;

for ifeature = 1:2
    xaxis = xaxis_tuning{ifeature};
    nfilters = length(xaxis);
    xaxis_itp = xaxis_itp_tuning{ifeature};
    ifamily = ifamily_perF(ifeature);
    nameParams = namesParams_all{ifamily};
    nparams_full = length(nameParams);
    namesParams_unit = namesParams_unit_perF{ifeature};
    
    for itype = 3
        for iparam = 1:nparams_full
            figure('Position', [0 200 ncomb6*300 500]);
            
            % get boostrapped data of all subj
            if nB==1, params = squeeze(margParams{ifeature}(:, itype, :, :, iparam)); % nsubj x nB x ncomb6 x nlines
            else, params = squeeze(margParams{ifeature}(:, :, itype, :, :, iparam)); % nsubj x nB x ncomb6 x nlines
            end
            
            % get median and SEM
            if nsubj == 1
%                 if nB==1
                    params_med_ = params;
                    SEM_neg = zeros(ncomb6, 2);
                    SEM_pos = SEM_neg;
                else, [params_med_, ~, ~, SEM_neg, SEM_pos] = getCI(params, 1, 1);
%                 end
                params_groupAVE = params_med_;
                params_groupSEM = zeros(size(params_groupAVE));
%             else % muted because the input margParams is already the median
%                 [params_med_, ~, ~, SEM_neg, SEM_pos] = getCI(params, 1, 2);
%                 params_groupAVE = squeeze(mean(params_med_));
%                 params_groupSEM = squeeze(std(params_med_)/sqrt(nsubj));
            end
            
            % and save the value for corrlation analysis
            params_all{ifeature, itype, iparam} = params;
            params_med_all{ifeature, itype, iparam} = params_med_;
            params_SEM_neg_all{ifeature, itype, iparam} = SEM_neg;
            params_SEM_pos_all{ifeature, itype, iparam} = SEM_pos;
            params_groupAVE_all{ifeature, itype, iparam} = params_groupAVE;
            params_groupSEM_all{ifeature, itype, iparam} = params_groupSEM;
            
            if ifeature == 1 % ORI
                ydiffErr = [.18, 16.5, 0.005];
                switch iparam
                    case 1, ticks_scatter = .05:.05:.25; ylim_scatter = [.04, .26]; ticks_bar = .09:.03:.21; ylim_bar = [.084, .216];
                    case 2, ticks_scatter = 10:2:18; ylim_scatter = [9.6, 18.6];ticks_bar = 10:2:18; ylim_bar = [9.6, 18.4];
                    case 3, ticks_scatter = -.04:.02:.04; ylim_scatter = [-.044, .044]; ticks_bar = -.02:.01:.02; ylim_bar = [-.022, .022];
                end
            else % SF
                ydiffErr = [2.05, 0.11 0.37, 0.01, 3.2];
                switch iparam
                    case 1, ticks_scatter = 1.4:.2: 2.2; ylim_scatter = [1.35, 2.25];
                    case 2, ticks_scatter = .03:.03:.15; ylim_scatter = [.024, .156]; ticks_bar = .04:.02:.12; ylim_bar = [.036, .124];
                    case 3, ticks_scatter = 1:.1:1.4; ylim_scatter =log2([.98 1.42]);
                    case 4, ticks_scatter = -.02:.01:.02; ylim_scatter = [-.022, .022];
                    case 5, ticks_scatter = 0:1.5:6; ylim_scatter = [-.7,  6.3];
                end
            end
            
            
            for icomb = 1:ncomb6
                
                if paramInd_perF{ifeature}(iparam) == 1, nlines = 2; else, nlines = 1; end
                
                %%%%%%%%
                %       bar       %
                %%%%%%%%
                subplot(2, ncomb6, icomb), hold on
                RCplot_param_bar(params_all, params_groupAVE_all, params_SEM_neg_all, ifeature, itype, iparam, icomb, iline, params_groupSEM, ticks_bar)
                
                %%%%%%%%
                %    scatter   %
                %%%%%%%%
                subplot(2, ncomb6, icomb+1), hold on, box on
                RCplot_params_scatter(ylim_scatter, params_med_, ifeature, ticks_scatter, iparam, params_SEM_neg, params_SEM_pos, params_groupAVE, params_groupSEM)
                
            end % end of icomb
            
            set(findall(gcf, '-property', 'FontSize'), 'FontSize',18)
            sgtitle(sprintf('[%s] %s %s', namesType{itype}, namesFeature{ifeature}, namesParams_unit{iparam}), 'FontSize',20)
            
        end % end of iparam
    end % end of itype
end % end of ifeature

%% helper function [bar]
function RCplot_param_bar(params_all, params_groupAVE_all, params_SEM_neg_all, ifeature, itype, iparam, icomb, iline, params_groupSEM, ticks_bar)
load('params_RCplot_temp.mat') % load some common parameters

sz_label = 12;
sz_ticks = 12;
wd_border = 1;
sz_marker = 10;

margParams_ave2 = nan(1,2);
for iline = 1:nlines
    
    color = colors_comb(combInd(icomb, iline), :);
    x = params_groupAVE_all{ifeature, itype, iparam}(icomb, iline);
    
    bar(iline, x, 'FaceColor', 'w', 'EdgeColor', color, 'BarWidth', .5, 'linewidth', wd_border)
    if nsubj==1
        error_neg = params_SEM_neg_all{ifeature, itype, iparam}(icomb, iline);
        error_pos = params_SEM_neg_all{ifeature, itype, iparam}(icomb, iline);
        errorbar(iline, x, error_neg, error_pos, '.', 'color', color, 'CapSize', 0, 'linewidth', wd_border)
    else
        error = params_groupSEM{ifeature}(iline, iparam);
        errorbar(iline, x, error, '.', 'color', color, 'CapSize', 0, 'linewidth', wd_border)
    end
end

yticklabels(ticks_bar)
if nlines == 2
    xlim([.3, 2.7])
    xticks(1:nlines)
    xticklabels({namesLocComb{combInd(icomb, :)}})
else
    xlim([.5, 1.5])
    xticks(1)
    xticklabels([namesLocComb{combInd(icomb, 1)}, ' and ', namesLocComb{combInd(icomb, 2)}])
end
ax = gca;
ax.XAxis.FontSize = sz_label;
ax.YAxis.FontSize = sz_ticks;
ax.LineWidth = wd_border;
box on

% plot extra lines
%                 switch ifeature
%                     case 1 % ORI (y=0 for baseline)
%                         if iparam == 3, yline(0, 'color', ones(1,3)*.5, 'linewidth', wd_border); end
%                     case 2 % SF (y=2 for peak; y=0 for baseline)
%                         if iparam == 1, yline(2, 'color', ones(1,3)*.5, 'linewidth', wd_border);
%                         elseif iparam == 4, yline(0, 'color', ones(1,3)*.5, 'linewidth', wd_border);
%                         end
%                 end
% ylabel
%                 ylabel(nameParams_unit{iparam}, 'fontsize', sz_label)

% compare two params
if nsubj==1 % compare across boots
    diff_sem = 0;
    if nB==1
        x1 = squeeze(params_all{ifeature, itype, iparam}(icomb, 1));
        x2 =  squeeze(params_all{ifeature, itype, iparam}( icomb, 2));
    else
        x1 = squeeze(params_all{ifeature, itype, iparam}(:, icomb, 1));
        x2 =  squeeze(params_all{ifeature, itype, iparam}(:, icomb, 2));
    end
else % compare across observers
    x1 = params_med{ifeature}(:, 1, iparam);
    x2 = params_med{ifeature}(:, 2, iparam);
    diff_sem = std(x1-x2)/sqrt(nsubj);
end

[~, p,~, ~] = ttest(x1, x2);

%                 errorbar(1.5, ydiffErr(iparam), diff_sem, 'k', 'CapSize', 0, 'linewidth', wd_border)
%                 plot([1,2], [ydiffErr(iparam), ydiffErr(iparam)], 'k-', 'linewidth', wd_border)
%                 string_s = getString_starts(p);
%                 text(1.5, ydiffErr(iparam)+diff_sem*1.5, string_s, 'HorizontalAlignment', 'center', 'fontsize', sz_title)

axis square
end

%% helper function [scatter]
function RCplot_params_scatter(ylim_scatter, params_med, ifeature, ticks_scatter, iparam, params_SEM_neg, params_SEM_pos, params_groupAVE, params_groupSEM)
load('params_RCplot_temp.mat') % load some common parameters

sz_marker = 20;
sz_label = 20;
sz_ticks = 15;
wd_border = 2;

% diagonal line
plot(ylim_scatter, ylim_scatter, 'k-', 'linewidth', wd_border)

% IDVD data with CI errorbar
x_F = params_med{ifeature}(:, 1, iparam);
x_P = params_med{ifeature}(:, 2, iparam);
x_P_SEMneg = params_SEM_neg{ifeature}(:, 2, iparam);
x_P_SEMpos = params_SEM_pos{ifeature}(:, 2, iparam);
x_F_SEMneg = params_SEM_neg{ifeature}(:, 1, iparam);
x_F_SEMpos = params_SEM_pos{ifeature}(:, 1, iparam);

errorbar(x_F, x_P, x_P_SEMneg, x_P_SEMpos, x_F_SEMneg, x_F_SEMpos, '.k', 'CapSize', 0, 'linewidth', wd_border/2)
%             errorbar(x_F, x_P, x_F_SEMneg, x_F_SEMpos, 'horizontal', '.r', 'CapSize', 0, 'linewidth', wd_border/2)
%             errorbar(x_F, x_P, x_P_SEMneg, x_P_SEMpos, 'vertical',  '.b', 'CapSize', 0, 'linewidth', wd_border/2)
scatter(x_F, x_P, sz_marker, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'linewidth', 2)

% group verage
x_F = params_groupAVE{ifeature}(1, iparam);
x_P = params_groupAVE{ifeature}(2, iparam);
x_F_SEM = params_groupSEM{ifeature}(1, iparam);
x_P_SEM = params_groupSEM{ifeature}(2, iparam);

errorbar(x_F, x_P, x_F_SEM, 'horizontal', '.','color', [.75,0,.75], 'CapSize', 0, 'linewidth', 8)
errorbar(x_F, x_P, x_P_SEM, 'vertical', '.','color', [.75,0,.75], 'CapSize', 0, 'linewidth', 8)

% limit
xlim(ylim_scatter)
ylim(ylim_scatter)

if (ifeature == 2) && (iparam == 3)
    xticks(log2(ticks_scatter))
    yticks(log2(ticks_scatter))
else
    xticks(ticks_scatter)
    yticks(ticks_scatter)
end
xticklabels(ticks_scatter)
yticklabels(ticks_scatter)

ax = gca;
ax.XAxis.FontSize = sz_ticks;
ax.YAxis.FontSize = sz_ticks;
ax.LineWidth = wd_border/1.2;

xlabel(namesLocComb{combInd(icomb, 1)}, 'FontSize', sz_label)
ylabel(namesLocComb{combInd(icomb, 2)}, 'FontSize', sz_label)

axis square

end

