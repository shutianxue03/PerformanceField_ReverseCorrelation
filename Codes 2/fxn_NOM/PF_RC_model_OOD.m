% estimate internal-external noise ratio by minimizing the calculated pA and estimated pA
% only analyze data, not plotting

clc
clear all
warning off
format compact

addpath(genpath('Data_OOD'))
addpath(genpath('fxn_model'))
addpath(genpath('fxn_analysis_RC'))

nsubj = 12; nsubj_real = nsubj;
iLoc_all = [1,8];
nLoc = length(iLoc_all);
nparams_model = 3; % 1 = alpha (assumed the same for PRS and ABS)
% 2 = alpha (assumed the same for PRS and ABS) and thresh
% 3 = alpha_PRS, alpha_ABS and thresh
fitNoiseMode.i = 2; % (1) metrics (2) trial-wise responses
ni = 100;
nrep = 20; % 20
namesMetrics = {'dprime',  'Criterion', 'pC', 'pHit', 'pFA', 'pA', 'pA [PRS]', 'pA [ABS]'};
nmetrics = length(namesMetrics);
markers_allSubj_full = {'o', 's', 'd', '^','v',  '<', '>','p', 'h', '+', 'x', '-'}; % for each subj
markers_allSubj = markers_allSubj_full(1:nsubj);

templateType.i = 3; % (1) raw kernel (2) positive kernel (3) reconstructed kernel
convolveType.i = 2; % (1) dot multiply (2) comvolution
IVType.i = 3; % (1) sum up (2) max
normIV = 0;

errorFxn = '@(pred, data) sum(-log(normpdf(pred, data, 1)), ''all'')';
errorComp = [1,2, 6:8]; % (1) dprime, (2) criterion, (3) pC, (4) pHit, (5) pFA, (6) pA, (7) pA_PRS, (8) pA_ABS

%% empty containers (_allSubj)
data_allSubj = nan(nsubj, nLoc, nmetrics);
pred_allSubj = nan(nsubj, nLoc, nmetrics, ni);
params_est_allSubj = nan(nsubj, nLoc, nparams_model, ni);

%% save the description of the analysis
templateType.all = {'(1) the 2D kernel (raw)',...
    '(2) the 2D kernel (all positive)', ...
    '(3) the reconstructed 2D kernel (all positive)', ...
    '(4) the energy profile of the gabor (CST=1)'};
convolveType.all = {'(1) template(:) .* energy(:)', '(2) conv2(template, energy)'};
IVType.all = {'(1) sum all channels up', '(2) select the highest value', '(3) normalization'};
fitNoiseMode.all = {'(1) metrics','(2) trial-wise responses'};

%% shell
for isubj = 1:nsubj % the selected subj among {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'AS', 'HL', 'FH', 'HA', 'CS', 'DT'};
    for iiLoc = 1:nLoc
        iLoc = iLoc_all(iiLoc);
        
        %% PARAMS
        subjList = {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'AS', 'HL', 'FH', 'HA', 'CS', 'DT'};
        nblocks_allSubj = [200, 240, 210, 240, 220, 220, 215, 220, 205, 210, 205, 205];
        
        subjName = subjList{isubj};
        nblocks = nblocks_allSubj(isubj);
        
        nB = 1;
        ntrials_sim = 1e5;
        nfiltersOri = 29;
        nfiltersSF = 29;
        itype = 3;
        
        nLoc = length(iLoc_all);
        options = optimoptions('fmincon','MaxIterations', 1e4, 'Display','off');
        fprintf('\nNumber of boostraps (nB): %d\nNumber of ORI filters: %d\nNumber of SF filters: %d\n\n', nB, nfiltersOri, nfiltersSF)
        
        fileName_idvd = sprintf('Data_OOD/%s_L%d_model.mat', subjName, iLoc);
        fileDir_idvd = dir(fileName_idvd);
        
        fprintf('\nSubj name: %s \nNumber of blocks (nblocks): %d\n', subjName, nblocks)
        
        
        %%  define the range
        % Params:
        % 1. internal-external noise ratio of the PRS trials
        % 2. internal-external noise ratio of the ABS trials
        % 3. criterion (not the criterion between -1 to 1)
        
        % Define the initial/lb/ub
        switch nparams_model
            case 1
                params0 = .5;
                params_lb = 0;
                params_ub = 2;
            case 2
                params0 = [.5, 0];
                params_lb = [0, -1];
                params_ub = [2, 1];
            case 3
                params0 = [.5, .5, .5];
                params_lb = [0, 0, 0];
                params_ub = [2, 2, 1];
        end
        assert(nparams_model == length(params_ub));
        
        load('Data_OOD/filtersFile', 'energy2D_gabor')
        
        %% create README
        README = sprintf('\n Loc fitted = %s\n ni = %d\n nrep = %d\n Template: %s\n Convolve type = %s\n Sum up or select the max = %s\n Normalize IV distributions (1=YES, 0=NO): %d\n Fit mode (1 = metrics, 2 = trial-wise resp): %d\n Criterion [initial/lb/ub]: %.1f, %.1f, %.1f\n Error fxn: %s\n Error components: %s\n \n', ...
            num2str(iLoc_all), nB, nrep, ...
            templateType.all{templateType.i}, ...
            convolveType.all{convolveType.i}, ...
            IVType.all{IVType.i}, ...
            normIV, ... % whether IV distributions are normalized?
            fitNoiseMode.all{fitNoiseMode.i}, ...
            [params0(1,3), params_lb(1,3), params_ub(1,3)], ...
            errorFxn,...
            num2str(errorComp)); % error components
        
        subjName_ = subjName;
        
        %%
        if isempty(fileDir_idvd)
            load(sprintf('Data_OOD/%s%d_B%d_%d_%d.mat', subjName_, nblocks, nB, nfiltersOri, nfiltersSF))
            load(sprintf('Data_OOD/%s%d_energy_%d_%d.mat', subjName_, nblocks, nfiltersOri, nfiltersSF))
            load(sprintf('Data_OOD/%s%d_behavMeas.mat', subjName_, nblocks))
            
            nAllTrials = length(energy2D_allT_perComb{1});
            ntrials = nAllTrials/2;
            
            subjName = subjName_;
            
            %% 1. calculate the internal variable
            fprintf('1. Calculating the internal variable... ')
            clear data % as the size of IV differs across subjects
            
            % get template and energy of that location 
            [template, e2D] = fxn_getTemplate(kernels2D_perComb, itype, iLoc, templateType, energy2D_allT_perComb);
            
            if IVType.i ~= 3
                %%%%%%%
                fxn_getIV
                %%%%%%%
                % plot IV distribution
                figure, hold on
                histogram(IV_PRS)
                histogram(IV_ABS)
                % plot the convoluted IV per channel
                data.imax_allT = imax_allT;
                figure % [imax_ORI_PRS, imax_ORI_ABS; imax_SF_PRS, imax_SF_ABS];
                subplot(2,2,1), histogram(squeeze(imax_allT(:, 1, 1))), title('[PRS] ORI'), xlim([1,size(template,1)]), xline(ORI_bound(1), 'r-', 'linewidth', 2); xline(ORI_bound(2), 'r-', 'linewidth', 2);
                subplot(2,2,2), histogram(squeeze(imax_allT(:, 1, 2))), title('[ABS] ORI'), xlim([1,size(template,1)]), xline(ORI_bound(1), 'r-', 'linewidth', 2); xline(ORI_bound(2), 'r-', 'linewidth', 2);
                subplot(2,2,3), histogram(squeeze(imax_allT(:, 2, 1))), title('[PRS] SF'), xlim([1,size(template,2)])
                subplot(2,2,4), histogram(squeeze(imax_allT(:, 2, 2))), title('[ABS] SF'), xlim([1,size(template,2)])
                sgtitle('At which ORI/SF channel max IV is found')
                % compile data
                data.IV(1, :) = IV_PRS; % before adding internal noise
                data.IV(2, :) = IV_ABS;
                
                fprintf('Done!\n')
            end
            
            % compile data
            data.resp(1,:) = data_both_allT_perComb{iLoc, 2}(1:ntrials); % #2 is the response
            data.resp(2,:) = data_both_allT_perComb{iLoc, 2}((ntrials+1):(ntrials*2)); % #2 is the response
            data.metrics_sim = [dprime_allT_perComb(iLoc), criterion_allT_perComb(iLoc), ...
                pC_allT_perComb(iLoc, [3,1,2]), pA_allT_perComb(iLoc, [3,1,2])];
            % keep the field name 'metrics_sim'!!
  
            
            
            %% 2. fit & predict
            fprintf('2. Fitting...')
            
            params_est_allB = nan(nparams_model, ni);
            pred_allB = nan(nmetrics, ni);
            
            tStart = tic;
            for ii = 1:ni
                
                if IVType.i == 3, fxn_getIV, end
                %%%%%%%
                %   Fitting  %
                %%%%%%%
                fxn_estParams = @(params) fxn_getError(fitNoiseMode.i, params, data, errorComp);
                problem_ML = createOptimProblem('fmincon','objective', fxn_estParams,'x0', params0, 'lb', params_lb, 'ub', params_ub, 'options',options);
                ms_ML = MultiStart('StartPointsToRun', 'bounds', 'Display', 'off', 'UseParallel', 1);
                
                params_est = run(ms_ML, problem_ML, nrep);
                %                 params_est = fmincon(fxn_estParams, params0, [], [], [], [], params_lb, params_ub, [], options);
                params_est_allB(:, ii) = params_est;
                
                %%%%%%%%%
                %   Prediction  %
                %%%%%%%%%
                pred = MR_pred(params_est, data);
                pred_allB(:, ii) = pred.metrics;
                
                if ~mod(ii, ni/10), fprintf('='), end
            end
            
            tEnd = toc(tStart);
            fprintf('Done! %.1f mins\n', tEnd/60)
            
            %  4. save IDVD data per loc
            save(fileName_idvd, 'template*', 'README', 'data','params_est_allB', 'pred_allB')
            fprintf(' ========== Finished ==========\n\n')
            
        else
            load(fileName_idvd)
        end
        data_allSubj(isubj, iiLoc, :) = data.metrics_sim;
        pred_allSubj(isubj, iiLoc, :, :) = pred_allB;
        params_est_allSubj(isubj, iiLoc, :, :) = params_est_allB;
        
        %% quick plot
        %
        %     figure('Position', [0 200 1500 200])
        %     for im = 1:nmetrics
        %         pred_im = squeeze(pred_allB(:, :, im));
        %         subplot(2,nmetrics, im), hold on
        %         histogram(pred_im)
        %         xline(data.metrics_sim(im), 'k-', 'linewidth', 2);
        %         title(namesMetrics{im})
        %
        %         subplot(2,nmetrics, im+nmetrics), hold on
        %         bar(data.metrics_sim(im))
        %         errorbar(median(pred_im), std(pred_im), 'o')
        %         if find(im == [3,4,6:8]), ylim([.5, 1])
        %         elseif im == 1,  ylim([0, 2])
        %         elseif im == 2, ylim([-1, 1])
        %         elseif im == 5,  ylim([0, 1])
        %         end
        %     end
        %     sgtitle(subjName_)
        %
        %     figure('Position', [0 200 600 200])
        %     for ip = 1:nparams_model
        %         p = squeeze(params_est_allB(:, :, ip));
        %         subplot(2,nparams_model, ip), hold on
        %         histogram(p)
        %
        %         subplot(2,nparams_model, ip+nparams_model), hold on
        %         errorbar(median(p), std(p), 'o')
        %     end
    end % end of iLoc
end % end of isubj

%% plot all subjects
nsubj = nsubj_real;
faceAlpha = 0.3; % for error band
marks_allSubj_full = {'o', 's', 'd', '^','v',  '<', '>','p', 'h', '+', 'x', '-'}; % for each subj
marks_allSubj = marks_allSubj_full(1:nsubj);
assert(length(marks_allSubj) == nsubj)

%% color
colorsType = [1,0,0;0,0,1;0,0,0];
colors_comb = [
    0,0,0; ...,     % center; black
    0, .75, 0; ..., % Left: light green
    1, 0, 0,; ...,   % upper: red
    0, .35, 0; ..., % right: dark green
    0, 0, 1; ...,    % lower: blue
    0, .5, 0; ...,   % HM: green
    .5, 0, 1; ...,    % VM: purple
    .5, .5, .5];      % peri: darker grey

colors5Loc = colors_comb(1:5, :);

%% plot measured vs.predicted metrics
figure('Position', [2000 200 1500 1500])
pred_allSubj_med = median(pred_allSubj, 4);
for im = 1:nmetrics
    pp = squeeze(pred_allSubj_med(:, :, im));
    dd = squeeze(data_allSubj(:, :, im));
    axis_max = max([pp(:); dd(:)]);
    axis_min = min([pp(:); dd(:)]);
    
    for iiLoc = 1:nLoc
        iLoc = iLoc_all(iiLoc);
        data_im = squeeze(data_allSubj(:, iiLoc, im));
        pred_im = squeeze(pred_allSubj_med(:, iiLoc, im));
        
        subplot(3,3, im), hold on
        for isubj = 1:nsubj
            plot(data_im(isubj), pred_im(isubj), markers_allSubj{isubj}, 'MarkerFaceColor', colors_comb(iLoc, :), 'MarkerEdgeColor', 'w')
        end
        errorbar(mean(data_im), mean(pred_im), std(data_im)/sqrt(nsubj-1), 'horizontal', 'linewidth', 2, 'color', colors_comb(iLoc, :), 'CapSize', 0)
        errorbar(mean(data_im), mean(pred_im), std(pred_im)/sqrt(nsubj-1), 'vertical', 'linewidth', 2, 'color', colors_comb(iLoc, :), 'CapSize', 0)
    end
    axis square, box on
    ylabel('Prediction')
    xlabel('Data')
    xlim([axis_min, axis_max])
    ylim([axis_min, axis_max])
    plot([axis_min, axis_max], [axis_min, axis_max], 'k-')
    title(namesMetrics{im})
    
end

%% plot estimated parameters (fovea vs. peri, PRS vs. ABS)
p1_PRS = squeeze(params_est_allSubj(:, 1, 1, :)); p1_PRS_med = median(p1_PRS, 2);
p1_ABS = squeeze(params_est_allSubj(:, 1, 2, :)); p1_ABS_med = median(p1_ABS, 2);
p2_PRS = squeeze(params_est_allSubj(:, 2, 1, :)); p2_PRS_med = median(p2_PRS, 2);
p2_ABS = squeeze(params_est_allSubj(:, 2, 2, :)); p2_ABS_med = median(p2_ABS, 2);

axis_min = min([p1_PRS_med;p1_ABS_med; p2_PRS_med;p2_ABS_med]);
axis_max = max([p1_PRS_med;p1_ABS_med; p2_PRS_med;p2_ABS_med]);

figure('Position', [2000 200 800 200])

for iplot = 1:2
    subplot(1,2,iplot), hold on, box on
    for ind = 1:2
        switch iplot
            case 1 % compare fovea (x) vs. peri (y)
                colors = colorsType;
                xlabel_ = 'Fovea';
                ylabel_ = 'Periphery';
                switch ind
                    case 1, x = p1_PRS_med;  y = p2_PRS_med;
                    case 2, x = p1_ABS_med; y = p2_ABS_med;
                end
            case 2 % compare PRS (x) vs. ABS (y)
                colors = colors_comb(iLoc_all, :);
                xlabel_ = 'PRS';
                ylabel_ = 'ABS';
                switch ind
                    case 1, x = p1_PRS_med;  y = p1_ABS_med;
                    case 2, x = p2_PRS_med; y = p2_ABS_med;
                end
        end
        for isubj = 1:nsubj
            plot(x(isubj), y(isubj), markers_allSubj{isubj}, 'color', colors(ind, :));
        end
        errorbar(mean(x), mean(y), std(x)/sqrt(nsubj-1), 'horizontal', 'linewidth', 2, 'color', colors(ind, :), 'CapSize', 0)
        errorbar(mean(x), mean(y), std(y)/sqrt(nsubj-1), 'vertical', 'linewidth', 2, 'color', colors(ind, :), 'CapSize', 0)
    end
    xlabel(xlabel_)
    ylabel(ylabel_)
    xlim([axis_min, axis_max])
    ylim([axis_min, axis_max])
    plot([axis_min, axis_max], [axis_min, axis_max], 'k-')
    % title(namesMetrics{im})
    axis square
end
% %% compare PRS vs. ABS
% subplot(1,2,2), hold on, box on
% for isubj = 1:nsubj
%     plot(p1_PRS_med(isubj), p1_ABS_med(isubj), markers_allSubj{isubj}, 'color', colors_comb(1, :));
%     plot(p2_PRS_med(isubj), p2_ABS_med(isubj), markers_allSubj{isubj}, 'color', colors_comb(8, :));
% end
% ylabel('PRS')
% xlabel('ABS')
% xlim([axis_min, axis_max])
% ylim([axis_min, axis_max])
% plot([axis_min, axis_max], [axis_min, axis_max], 'k-')
% axis square
