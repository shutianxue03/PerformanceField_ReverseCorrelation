%% check the correlation between energy bins
% each bin is the energy component of a specific ORI and SF channel of all
% trials

load('Data_OOD/params_RC.mat')


%% SX's data
load('Data_OOD/AS215/AS_energy_29_29.mat')
load('Data_OOD/AS215/AS_behavMeas.mat')
iLoc = 2; % index of location (1=Fovea; 2=Left; 3=UVM; 4=Right; 5=LVM)
indLoc = dataMatrix(:, 5) == iLoc;
indiPRS = dataMatrix(:, 6); % 1= signal-PRS; 0=ABS

e3D = squeeze(e2D_allT(indLoc & ~indiPRS, :, :)); 
clear std
[ntrials, nORI, nSF] = size(e3D);

%% AF's data
% load('AF_EndoRC/source_energy.mat')
e3D_AF = patch(trialsMat(:, 3)==1 & trialsMat(:, 5)==1 & trialsMat(:, 9)==.26);
ntrials = length(e3D_AF);
[nORI, nSF] = size(e3D_AF{1});
e3D = nan(ntrials, nORI, nSF);
for itrial = 1:ntrials
    e3D(itrial,:,:) = e3D_AF{itrial};
end
axis_tuning = {1:nORI, 1:nSF};
axisTL_tuning = {1:5, 1:5};
axisTicks_tuning ={1:5, 1:5};

%%

assert(nORI<50 && nSF < 50); 

%% check 2D correlation between every pair of energy images
% takes a lot time!! to be more effcient, increase the interval of itrials_short
% itrials_short = 1:ntrials;
% ntrials_short = length(itrials_short);
% R = nan(ntrials_short, ntrials);
% for itrial = 1:ntrials_short
%     d1 = squeeze(e3D(itrial, :, :));
%     parfor itrial_c = 1:ntrials
%         d2 = squeeze(e3D(itrial_c, :, :));
%         R(itrial, itrial_c) = corr2(d1, d2);
%     end
% end
% figure, imagesc(R), colorbar


%% the main loop
close all
clc

params_est_SF_all = nan(nORI, nSF, 3);
params_est_ORI_all  = nan(nORI, nSF, 3);
Rsquared_SF_all = nan(nORI, nSF);
Rsquared_ORI_all = Rsquared_SF_all;
 
plotFlag = 1;

for iORI = 10%1:nORI
    for iSF = 10%1:nSF
        ebin_fix = squeeze(e3D(:, iORI, iSF));
        close all
        if plotFlag, figure('Position', [3e3 200 800 600]), end
        
        r_fix = nan(nORI, nSF);
        p_fix = r_fix;
        dotP_fix = r_fix;
        
        % standardization
        ebin_fix_stand = (ebin_fix - mean(ebin_fix))/std(ebin_fix);
        
        %
        for iORI_comp = 1:nORI
            for iSF_comp = 1:nSF
                ebin_comp = squeeze(e3D(:, iORI_comp, iSF_comp));
                % standardization
                ebin_comp_stand = (ebin_comp - mean(ebin_comp))/std(ebin_comp);
                
                % correlation 
                [r, p] = corr(ebin_fix_stand, ebin_comp_stand);
                r_fix(iORI_comp, iSF_comp) = r;
                p_fix(iORI_comp, iSF_comp) = p;
                
                % dot product (no standardization)
                dotP_fix(iORI_comp, iSF_comp) = sum(ebin_fix(:) .* ebin_comp(:));
                
            end % iSF_comp
        end % iORI_comp
        
        % plot SF and ORI marg
        % 1. Pearson's r
        [params_est_SF, params_est_ORI, Rsquared_SF, Rsquared_ORI] = fxn_plot([1,4,7], r_fix, iORI, iSF, 'Pearson''s r', axis_tuning, axisTL_tuning, axisTicks_tuning, plotFlag);
        if plotFlag
            % 2. p value
            fxn_plot([2,5,8], p_fix, iORI, iSF, 'P-value', axis_tuning, axisTL_tuning, axisTicks_tuning, plotFlag);
            % 3. dot product
            fxn_plot([3,6,9], dotP_fix, iORI, iSF, 'Dot product', axis_tuning, axisTL_tuning, axisTicks_tuning, plotFlag);
        end
        
        if plotFlag
            set(findall(gcf, '-property', 'FontSize'), 'FontSize',13)
            sgtitle(sprintf('Correlation between energy bins\nFixed bin: ORI = %d deg, SF = %.2f cpd', axis_tuning{1}(iORI), 2^axis_tuning{2}(iSF)))
        end
        if (Rsquared_SF<.5) || (Rsquared_ORI<.5)
            fprintf('iORI = %d, iSF = %d, R2_ORI = %.2f, R2_SF = %.2f\n', iORI, iSF, Rsquared_ORI, Rsquared_SF)
        end
        params_est_SF_all(iORI, iSF, :) = params_est_SF;
        params_est_ORI_all(iORI, iSF, :) = params_est_ORI;
        Rsquared_SF_all(iORI, iSF) = Rsquared_SF;
        Rsquared_ORI_all(iORI, iSF) = Rsquared_ORI;
        
        if plotFlag, waitforbuttonpress, end
        fprintf('=')
    end
        fprintf(' iORI = %d/%d  done\n', iORI, nORI)
end

fprintf('done\n')

%% plot the estimated parameters
namePs = {'Width', 'Gain', 'Baseline'};
np = 1;
for ip = 1:np
    subplot(2,np+1,ip),hold on
    slice_ORI = squeeze(params_est_ORI_all(:, :, ip));
    imagesc(axis_tuning{2}, axis_tuning{1}, slice_ORI), colorbar
    xticks(axisTicks_tuning{2}), xticklabels(axisTicks_tuning{2})
    yticks(axisTicks_tuning{1}), yticklabels(axisTicks_tuning{1})
    title(sprintf('ORI %s\n mean = %.2f, SD = %.2f', namePs{ip}, mean(slice_ORI(:)), std(slice_ORI(:))))
    
    subplot(2,np+1,ip+np+1),hold on
    slice_SF = squeeze(params_est_SF_all(:, :, ip));
    imagesc(axis_tuning{2}, axis_tuning{1}, slice_SF), colorbar
    xticks(axisTicks_tuning{2}), xticklabels(axisTicks_tuning{2})
    yticks(axisTicks_tuning{1}), yticklabels(axisTicks_tuning{1})
    title(sprintf('SF %s (octave) \nmean = %.2f, SD = %.2f', namePs{ip}, mean(slice_SF(:)), std(slice_SF(:))))
end

subplot(2,np+1,np+1), hold on
imagesc(Rsquared_ORI_all), colorbar
% xticks(axisTL_tuning{2}), xticks(axisTicks_tuning{2})

subplot(2,np+1,(np+1)*2), hold on
imagesc(Rsquared_SF_all), colorbar


%% helper function - plot
function [params_est_SF, params_est_ORI, Rsquared_SF, Rsquared_ORI] = fxn_plot(iplots, data, iORI, iSF, title_, axis_tuning, axisTL_tuning, axisTicks_tuning, plotFlag)

params_est_SF = nan;
params_est_ORI = nan;

%%%%%%%%%%%%
% 1. UPPER: colormap %
%%%%%%%%%%%%
if plotFlag
    subplot(3,3,iplots(1)), 
    imagesc(axis_tuning{2}, axis_tuning{1}, data), colorbar, xlabel('SF'), ylabel('ORI')
    % axis square
    xline(axis_tuning{2}(iSF), 'r-', 'linewidth', 2);
    yline(axis_tuning{1}(iORI), 'r-', 'linewidth', 2);
    yticks(axisTicks_tuning{1}), yticklabels(axisTL_tuning{1}), 
    xticks(axisTicks_tuning{2}), xticklabels(axisTL_tuning{2})
    title(title_)
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % 2. MIDDLE: marginalize across ORI (x is SF)
    %%%%%%%%%%%%%%%%%%%%%%%
    subplot(3,3, iplots(2)), hold on
    plot(axis_tuning{2}, data(iORI, :), 'ko'), xlabel('SF'), ylabel(title_)
    if ~ iplots(1) == 3, yticks(0:.2:1), yticklabels(0:.2:1), ylim([0,1]), end
    xline(axis_tuning{2}(iSF), 'r-', 'linewidth', 2);
    xticks(axisTicks_tuning{2}), xticklabels(axisTL_tuning{2})
end

% fit gaussian
if iplots(1) == 1
    params0 = [1,1,0]; params_lb = [.001, 0,  -1]; params_ub = [1, 2, .2];
    %     params0 = 1; params_lb = .001; params_ub = 1; % just fitting width, assuming gain=1 and base=0
    [params_est_SF, Rsquared_SF] = fxn_fitGaussian(data(iORI, :), axis_tuning{2}, iSF, params0, params_lb, params_ub, plotFlag);
end

%%%%%%%%%%%%%%%%%%%%%%%
% 3. LOWER: marginalize across SF (x is ORI)
%%%%%%%%%%%%%%%%%%%%%%%
if plotFlag
    subplot(3,3, iplots(3)), hold on
    plot(axis_tuning{1}, data(:, iSF), 'ko'), xlabel('ORI'), ylabel(title_)
    if ~ iplots(1) == 3, yticks(0:.2:1), yticklabels(0:.2:1), ylim([0,1]), end
    xline(axis_tuning{1}(iORI), 'r-', 'linewidth', 2);
    xticks(axisTicks_tuning{1}), xticklabels(axisTL_tuning{1})
end
% fit gaussian
if iplots(1) == 1
    params0 = [20, 1, 0]; params_lb = [.001, 0, -1]; params_ub = [45, 2, .2];
    %     params0 = 20; params_lb = .001; params_ub = 45;
    [params_est_ORI, Rsquared_ORI] = fxn_fitGaussian(data(:, iSF), axis_tuning{1}, iORI, params0, params_lb, params_ub, plotFlag);
end

end

%% helper fxn - predict gaussian
function y = predGaussian(x, params, center)
SX_normPDF = @(x,mu,sigma) exp(-(x-mu).^2/(2*sigma^2));

width = params(1);
gain = params(2);
base = params(3);
y = gain * SX_normPDF(x, center, width) + base;
end

% function y = predGaussian(x, width, center)
% SX_normPDF = @(x,mu,sigma) exp(-(x-mu).^2/(2*sigma^2));
% y = SX_normPDF(x, center, width);
% end

%% helper fxn - fit gaussian
function [params_est, Rsquared] = fxn_fitGaussian(y, x, ii, params0, params_lb, params_ub, plotFlag)
if size(y,1)>1, y = y'; end % make sure data has the same shape as x
options = optimoptions('fmincon','MaxIterations', 5000,'Display','off');
if params_ub(1) > 5 % x is ORI, to avoid the asymmetry playing confounding the fitting
    buffer = floor(length(y)/2)-2;
    if ii - buffer < 1, ii_start = 1; ii_end = ii+buffer;
    elseif  ii+buffer > length(y), ii_start = ii-buffer; ii_end = length(y);
    else, ii_start = ii-buffer; ii_end = ii+buffer;
    end
else, ii_start = 1; ii_end = length(y);
end

if plotFlag
    plot(x(ii_start:ii_end), y(ii_start:ii_end), 'k+') % plot the dots involved in fitting gaussian
end
center = x(ii);
fxn_getError_ = @(params) sumsqr(y(ii_start:ii_end) - predGaussian(x(ii_start:ii_end), params, center));
[params_est, error] = fmincon(fxn_getError_, params0, [], [], [], [], params_lb, params_ub, [], options);
ypred = predGaussian(x(ii_start:ii_end), params_est, center);
Rsquared = 1-sumsqr(y(ii_start:ii_end)-ypred)/sumsqr(y(ii_start:ii_end)-mean(y(ii_start:ii_end)));
if plotFlag
    plot(x(ii_start:ii_end), ypred,'c-')
    title(sprintf('gain = %.2f / width = %.2f / base = %.2f \nR^2 = %.2f', params_est, Rsquared))
end
end