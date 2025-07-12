clc
line_extra = [0,1]; % x=0 deg // x=1 cpd // x=1 cpd (zoomed in, delete this if unneeded)

combInd_inv = [1,1;4,1;3,2;4,2;3,1;2,1;2,2;1,2];

xaxis_itp_all = cell(1,2);

yticks_ = {-.05:.1:.25, -.02:.06:.16};
ylimit = {[-.07, .27], [-.032, .172]};
titleF = {'Ori', 'SF'};

%%
if nsubj == 1
    % when nsubj==1, the dot is median, the errorbar is 68% CI
    % otherwse, the dot is group ave, the errorbar is SEM
    marg_aveSubj_all = marg_medB_allSubj_all;
    marg_sem1Subj_all = marg_CI_B_pos_allSubj_all; 
    marg_sem2Subj_all = marg_CI_B_neg_allSubj_all; 
    pred_aveSubj_all = pred_medB_allSubj_all; 
    pred_sem1Subj_all = pred_CI_B_pos_allSubj_all;
    pred_sem2Subj_all = pred_CI_B_neg_allSubj_all;
else
    marg_sem1Subj_all = marg_semSubj_all;
    marg_sem2Subj_all = marg_semSubj_all;
    pred_sem1Subj_all = pred_semSubj_all;
    pred_sem2Subj_all = pred_semSubj_all;
end

%%
figure('Position', [0 600 1500 length(line_extra)*250])

for ifeature = 1:2 % 1=ori, 2=SF
    xaxis = axis_tuning{ifeature};
    nfilters = length(xaxis);
    xaxis_itp = axis_itp_tuning{ifeature};
    
    for iLoc = 1:nLoc5
        for itype = 1:ntypes
            
            % extract data
            marg5_aveSubj = marg_aveSubj_all{ifeature, itype, iLoc}.'; % 
            marg5_sem1Subj = marg_sem1Subj_all{ifeature, itype, iLoc}.';
            marg5_sem2Subj = marg_sem2Subj_all{ifeature, itype, iLoc}.';
            margPred5_aveSubj = pred_aveSubj_all{ifeature, itype, combInd_inv(iLoc, 1)}(combInd_inv(iLoc, 2), :);
            margPred5_sem1Subj = pred_sem1Subj_all{ifeature, itype, combInd_inv(iLoc, 1)}(combInd_inv(iLoc, 2), :);
            margPred5_sem2Subj = pred_sem2Subj_all{ifeature, itype, combInd_inv(iLoc, 1)}(combInd_inv(iLoc, 2), :);
            
            subplot(length(line_extra), nLoc5, nLoc5*(ifeature-1) + iLoc), hold on, box on
            
            % extra lines
            yline(0, 'handlevisibility', 'off', 'linewidth', 1, 'color', [.7, .7, .7]);
            xline(line_extra(ifeature), 'handlevisibility', 'off', 'linewidth', 1, 'color', [.7, .7, .7]);
            
            % plot raw data
            plot(xaxis, marg5_aveSubj, ['o', colorsType{itype}])
            errorbar(xaxis, marg5_aveSubj, marg5_sem1Subj, marg5_sem2Subj, ['.', colorsType{itype}], 'CapSize',0)
            
            % plot fitting
            plot(xaxis_itp, margPred5_aveSubj, ['-', colorsType{itype}])
            patch([xaxis_itp, flip(xaxis_itp)], [margPred5_aveSubj-margPred5_sem2Subj, flip(margPred5_aveSubj+margPred5_sem1Subj)], colorsType{itype}, 'FaceAlpha', .3, 'linestyle', 'none')
            
            % plot tuning curve of gabor and calculate the corr coeff
            %             marg5_gabor = reshape(marg5_gabor, length(marg5_gabor), 1);
            %             %             plot(xaxis, marg5_gabor, 'g-')
            %             [r, p] = corr([marg5_gabor, marg5_aveSubj']);
            %             string_stars = getString_starts( p(1,2));
            %             text_corr = sprintf('r=%.2f%s',r(1,2), string_stars);
            %             text(xticks_{ifeature}(1), y_max-.02*itype, text_corr, 'color', colorsType{itype})
            
            % labels and title
            if iLoc == 1, ylabel(titleF{ifeature}), xlabel(namesFeature{ifeature}), end
            title(namesLoc2D{iLoc})
            
            % ticks and limits
            xticks(axisTicks_tuning{ifeature})
            xticklabels(axisTL_tuning{ifeature})
            xlim(axisLim{ifeature})
            yticks(yticks_{ifeature})
            yticklabels(yticks_{ifeature})
            ylim(ylimit{ifeature})
            
        end % end of itype
    end % end of iLoc
end % end of ifeature

set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)
set(findall(gcf, '-property', 'linewidth'), 'linewidth',1.5)

make_sgtitle(sprintf('[%s] Marginalized values ', titleF{ifeature}), subjName, nB, nsubj, nAllTrials)

