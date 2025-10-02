
%%%%%%%%%%%%%%%% PART BEHAV %%%%%%%%%%%%%%%%
% focus on analysis of non-RC performance

close all, clc

addpath(genpath('Data_compile'))
addpath(genpath('fxn_analysis_RC_v2'))
addpath(genpath('XueCarrasco_JN/code'))

SX_RC1_setting

% directory
nameFigFolder_behav = 'XueCarrasco_JN/fig/BEHAV';
itype = 2;
nameFolderCompile_behav = 'Data_compile/BEHAV';

sz_h = 250; % height of the figure
sz_h_comp = sz_h/1.5; % height of the figure, compressed
load(sprintf('%s/n%d_B%d_perf_%s.mat', nameFolderCompile_behav, nsubj, nB, nameFileBehav))

%% 1. Plot perf of 3 locations
flag_plotIDVD = 0;
close all
sz_wd_perBar_comp = 60;
% iLocComb_all = [6,5,3];
% nBars = length(iLocComb_all);

for iPerf = 4%3:6
    switch iPerf
        case 1
            perf_allSubj = cs_allSubj;
            namePerf = 'CS';
            ticks_scatter = 1.5:.6:3.9;
            yline_ = nan;
            sz_fig = [nBars*sz_wd_perBar sz_h];
        case 2
            perf_allSubj = pA_allSubj;
            namePerf = 'pA';
            ticks_scatter = linspace(.6, .8, 5);
            yline_ = nan;
            sz_fig = [nBars*sz_wd_perBar sz_h];
        case 3
            perf_allSubj = dprime_allSubj;
            namePerf = 'dprime';
            ticks_scatter = 0:2;
            yline_ = dprime_theo;
            sz_fig = [nBars*sz_wd_perBar_comp sz_h_comp];
        case 4
            perf_allSubj = criterion_allSubj;
            namePerf = 'criterion';
            ticks_scatter = -1:1;
            yline_ = 0;
            sz_fig = [nBars*sz_wd_perBar_comp sz_h_comp];
        case 5
            perf_allSubj = RT_allSubj*1e3;% convert to ms
            namePerf = 'RT';
            ticks_scatter = 0:150:300;
            yline_ = 0;
            sz_fig = [nBars*sz_wd_perBar_comp sz_h_comp];
        case 6
            perf_allSubj = pC_allSubj*1e2;% convert to %
            namePerf = 'pC';
            ticks_scatter = 50:20:90;
            yline_ = 70;
            sz_fig = [nBars*sz_wd_perBar_comp sz_h_comp];
    end
    basicFxn_drawBars(squeeze(median(perf_allSubj, 2)), yline_, colors_comb(iLocComb_all, :), namesLocComb(iLocComb_all), ticks_scatter, round(ticks_scatter, 2), ...
        flag_plotIDVD, 0, namePerf, 0, sz_fig)
    saveas(gcf, sprintf('%s/n%d_%s.jpg', nameFigFolder_behav, nsubj, namePerf))
end

%% 2. plot HM vs. VM and L vs. UVM (CS)
sz_wd_perBar = 90;
nBars = 2;
sz_fig = [nBars*sz_wd_perBar, 300];

for ii=2
    switch ii
        case 1, data = cs_allSubj; xticks_ = linspace(1.6, 3.6, 5); dataName = 'CS';
        case 2, data = pA_allSubj; xticks_ = round(linspace(.5, .9, 5), 2); dataName = 'pA';
    end
    % data = dprime_allSubj; xticks_ = linspace(0,2, 5); dataName = 'dprime';
    
    HM = squeeze(data(:, :, 1));
    VM = squeeze(data(:, :, 2)+data(:, :, 3))/2;
    LVM = squeeze(data(:, :, 2));
    UVM = squeeze(data(:, :, 3));
    
    if dataName(1)=='N'
        HM = squeeze(data(:, 1, :));
        VM = squeeze(data(:, 2, :)+data(:, 3, :))/2;
        LVM = squeeze(data(:, 2, :));
        UVM = squeeze(data(:, 3, :));
    end
    
    basicFxn_drawBars([getCI(HM, 1, 2), getCI(VM, 1, 2)], [], colors_comb([6,7], :), {'HM', 'VM'}, xticks_, xticks_, flag_plotIDVD, 1, 'Asymmetry', 0, sz_fig);
    saveas(gcf, sprintf('%s/n%d_%s_HM_VM.jpg', nameFigFolder_behav, nsubj, dataName))
    
    basicFxn_drawBars([getCI(LVM, 1, 2), getCI(UVM, 1, 2)], [], colors_comb([5,3], :), {'LVM', 'UVM'}, xticks_, xticks_, flag_plotIDVD, 1, 'Asymmetry', 0, sz_fig);
    saveas(gcf, sprintf('%s/n%d_%s_LVM_UVM.jpg', nameFigFolder_behav, nsubj, dataName))
end

%% for FOV vs PERI
% xticks_ =.6:.05:.8; dataName = 'pA';
% basicFxn_drawBars(mm_med, [], colors_comb([1,8], :), {'Fov', 'Peri'}, xticks_, xticks_, flag_plotIDVD, 1, 'Asymmetry', 0, sz_fig);
% saveas(gcf, sprintf('%s/n%d_%s_L18.jpg', nameFigFolder_behav, nsubj, dataName))
%

%% 3. plot HVA vs.VMA
% close all
nBars = 2;
sz_fig = [nBars*sz_wd_perBar_comp, 200];
HVA = (HM-VM)./(HM+VM)*100;
VMA = (LVM-UVM)./(LVM+UVM)*100;
xticks_= -5:5:15;
basicFxn_drawBars([getCI(HVA, 1, 2), getCI(VMA, 1, 2)], [], [1,1,1;1,1,1], {'HVA', 'VMA'}, xticks_, xticks_, 1,1, 'Asymmetry', 0, sz_fig);
% saveas(gcf, sprintf('%s/n%d_%s.jpg', nameFigFolder_behav, nsubj, 'CS_asym'))

