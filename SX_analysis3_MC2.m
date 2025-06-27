close all
clc
addpath(genpath('Data_MC2'))

SX_RC1_setting
MCmode = 3;

nORI=37;
nSF=37;
text_m='m';
min_f_perF = cell(2,1);
ifamily_all_perF = {[1,8], [2,3,12]};
iLocComb_all = 2:7;
figFolderName = '';
nIC=1; if MCmode==3, nIC=3; end
nblocks_allSubj = [200, 240, 210, 240, 220, 220, 215, 220, 205, 210, 205, 205];
colors = {'c', 'm', 'g'};
subjList = {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'AS', 'HL', 'FH', 'HA', 'CS', 'DT'};
nblocks_allSubj = [200, 240, 210, 240, 220, 220, 215, 220, 205, 210, 205, 205];

indSubj = [1:10, 12];
subjList = subjList(indSubj);
nblocks_allSubj = nblocks_allSubj(indSubj);
markers_allSubj = markers_allSubj(indSubj);
nsubj = length(indSubj);
% min_f_allSubj = nan(nsubj, 2, 8, nIC);
if MCmode<3, dev_allSubj = nan(nsubj, 12, 8); else, dev_allSubj = nan(nsubj, 12, 8, nIC); end

for isubj = 1:nsubj
    subjName = subjList{isubj};
    
    % load MC data
    load(sprintf('Data_MC2/%s_%d%d%s_mode%d', subjName, nORI, nSF, text_m, MCmode))
    if MCmode<3
        dev_allSubj(isubj, :, :) = dev_allF_allLoc(:, :, 1); % onlylook at median
    else
        dev_allSubj(isubj, :, :, :) = dev_allF_allLoc;
    end
    
    for ifeature=1:2
        ifamily_all = ifamily_all_perF{ifeature};
        figFolderName = sprintf('XueCarrasco/fig/MC/singleLoc/%d%d%s/mode%d/%s',nORI, nSF, text_m, MCmode, namesFeature{ifeature});
        folderDir = dir(figFolderName); if isempty(folderDir), mkdir(figFolderName), end
        
        for iIC=1:nIC
            dd = squeeze(dev_allF_allLoc(ifamily_all, :, iIC)-min(dev_allF_allLoc(ifamily_all, :, iIC)));
            for iLocComb = iLocComb_all
                figure('Position', [0 0 1e3 400])
                % dev of all families
                subplot(1, 2, 1), hold on
                for ifamily = ifamily_all
                    bar(find(ifamily == ifamily_all), dd(ifamily == ifamily_all, iLocComb), ...
                        'EdgeColor', colors{ifamily==ifamily_all}, 'FaceColor', colors_comb(iLocComb, :))
                end
                [~, min_f] = min(dd(:, iLocComb));
                min_f = ifamily_all(min_f);
                %                 min_f_allSubj(isubj, ifeature, iLocComb, iIC) = min_f;
                title(sprintf('F%d wins', min_f))
                xticks(1:length(ifamily_all))
                xticklabels(namesFamily_all(ifamily_all))
                % data and pred of all families
                subplot(1, 2, 2), hold on
                plot(axis_tuning{ifeature}, yData_allLoc{iLocComb, ifeature}, 'o', 'color', colors_comb(iLocComb, :), 'handlevisibility', 'off')
                for ifamily = ifamily_all
                    plot(axis_tuning{ifeature}, margPred_allF_allLoc{ifamily, iLocComb}, ['-', colors{ifamily==ifamily_all}])
                end
                legend(namesFamily_all(ifamily_all), 'Location', 'best')
                
                set(findall(gcf, '-property', 'LineWidth'), 'LineWidth',2)
                set(findall(gcf, '-property', 'fontsize'), 'fontsize',20)
                
                if MCmode==3
                    sgtitle(sprintf('%s %s [%s] %s', subjName, namesLocComb{iLocComb}, namesFeature{ifeature}, namesIC{iIC}))
                    saveas(gcf, sprintf('%s/%s_L%d_IC%d.jpg', figFolderName, subjName, iLocComb, iIC))
                else
                    sgtitle(sprintf('%s %s [%s]', subjName, namesLocComb{iLocComb}, namesFeature{ifeature}))
                    saveas(gcf, sprintf('%s/%s_L%d.jpg', figFolderName, subjName, iLocComb))
                end
            end % iLocComb
        end % iIC
        close all
    end % ifeature
    
end % isubj


%%
iLocComb_all=[6,5,3,7];
iIC = 1;

for ifeature=2
    if ifeature==1,ifamily_all = [1,8];
    else, ifamily_all = [2,3,12];
    end
    if MCmode<3
        dev_allSubj_ = dev_allSubj(:, ifamily_all, iLocComb_all);
    else
        dev_allSubj_ = dev_allSubj(:, ifamily_all, iLocComb_all, iIC);
    end
    dev_allSubj_delta = dev_allSubj_-min(dev_allSubj_(:));
    [dev_delta_ave, ~, ~, dev_delta_sem] = getCI(dev_allSubj_delta, 2, 1);
    
    figure('Position', [0 0 1.5e3 400])
    for iLocComb = iLocComb_all
        subplot(1, length(iLocComb_all), find(iLocComb==iLocComb_all)), hold on
        bar(1:length(ifamily_all), squeeze(dev_delta_ave(:, iLocComb==iLocComb_all)), 'FaceColor', 'w')
        for isubj = 1:nsubj
            plot(1:length(ifamily_all), dev_allSubj_delta(isubj, :,  iLocComb==iLocComb_all), ['-', markers_allSubj{isubj}], 'color', ones(1,3)/2)
        end
        xticks(1:length(ifamily_all))
        xticklabels(namesFamily_all(ifamily_all))
        title(namesLocComb{iLocComb_all(iLocComb==iLocComb_all)})
    end
    
    IV1 = nan(nsubj, length(ifamily_all), length(iLocComb_all));
    IV2=IV1;
    for isubj = 1:nsubj
        IV1(isubj, :, :) = repmat((1:length(ifamily_all))', 1,length(iLocComb_all));
        IV2(isubj, :, :) = repmat(1:length(iLocComb_all), length(ifamily_all), 1);
    end
    text_ANOVA = print_nANOVA({'Family', 'Loc'}, dev_allSubj_delta(:), {IV1(:), IV2(:)}, nsubj);
    if MCmode<3
        sgtitle(sprintf('n%d %s mode%d\n%s', nsubj, namesFeature{ifeature}, MCmode, text_ANOVA))
    else
        sgtitle(sprintf('n%d %s IC%d\n%s', nsubj, namesFeature{ifeature}, iIC, text_ANOVA))
    end
end