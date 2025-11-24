function SX_plot_NOM_kernels(nIterations, iModelA, iLocComb_all, subjList)
% SECTION 1: Template (kernel) reconstruction quality

SX_RC1_setting;

for s = 1:numel(subjList)
    subjName = subjList{s};
    for iLocComb = iLocComb_all
        % Folder and file naming must match compIV stage
        nameFolder_NOM = sprintf('%s/%s/L%d', nameFolder_Data_NOM_Trialwise, subjName, iLocComb);
        nameFile_compIV = sprintf('%s/n%d_A%d_compIV.mat', nameFolder_NOM, nIterations, iModelA);
        if ~isfile(nameFile_compIV)
            fprintf('Missing compIV file: %s\n', [nameFile_compIV '.mat']);
            continue;
        end

        % load recovered template and the true template
        load(nameFile_compIV, 'data_allB', 'template_true');

        % Extract recovered kernels and calculate correlation with ideal
        % template
        template_allIter = data_allB.template_allIter; % [nIterations x nORI x nSF]
        corr_iter = nan(nIterations,1);
        for ii = 1:nIterations
            kVec = squeeze(kernel2D_allB(ii,:,:));
            kVec = kVec(:);
            C = corrcoef(ideal_vec, kVec);
            corr_iter(ii) = C(1,2);
        end
        

        % 1.1 Subject × location mean kernel
        meanKernel = squeeze(mean(kernel2D_allB,1)); % [nORI x nSF]

        figure('Name', sprintf('Kernel_%s_L%d', subjName, iLocComb));
        imagesc(meanKernel); axis image xy;
        colorbar;
        title(sprintf('Reconstructed kernel: %s, %s', subjName, namesLocComb{iLocComb}));
        xlabel('SF'); ylabel('ORI');
        colormap(parula);

        % 1.2 Group-averaged per loc will be done after loop
        % (simplest: accumulate, then plot once per loc)
        % Here: store per-subject mean kernel
        outKernel.(sprintf('S%d_L%d', s, iLocComb)) = meanKernel; %#ok<STRNU>
        
        % 1.3 Correlation with ideal template per iteration
        ideal_vec = template_true(:);
        corr_iter = nan(nIterations,1);
        for ii = 1:nIterations
            kVec = squeeze(kernel2D_allB(ii,:,:));
            kVec = kVec(:);
            C = corrcoef(ideal_vec, kVec);
            corr_iter(ii) = C(1,2);
        end
        [corr_pt, corr_lo, corr_hi] = getCI(corr_iter); %#ok<*NASGU>

        % Optional: store for later group-level plot
        corr_store(s, iLocComb) = corr_pt; %#ok<AGROW>
    end
end

% 1.2 Group-averaged kernel per location (across subjects)
for iLocComb = iLocComb_all
    kernels_subj = [];
    for s = 1:numel(subjList)
        key = sprintf('S%d_L%d', s, iLocComb);
        if isfield(outKernel, key)
            kernels_subj = cat(3, kernels_subj, outKernel.(key));
        end
    end
    if isempty(kernels_subj), continue; end

    meanGroupKernel = mean(kernels_subj, 3);

    figure('Name', sprintf('GroupKernel_L%d', iLocComb));
    imagesc(meanGroupKernel); axis image xy;
    colorbar;
    title(sprintf('Group kernel, %s', namesLocComb{iLocComb}));
    xlabel('SF'); ylabel('ORI');
    colormap(parula);
end

% 1.3 Correlation summary figure
figure('Name','TemplateRecoveryCorr');
hold on;
for iLocComb = iLocComb_all
    vals = corr_store(:, iLocComb);
    vals = vals(~isnan(vals));
    if isempty(vals), continue; end
    mu = mean(vals);
    se = std(vals)/sqrt(numel(vals));
    x = find(iLocComb_all==iLocComb);
    errorbar(x, mu, se, 'o', 'MarkerFaceColor','k');
end
set(gca, 'XTick', 1:numel(iLocComb_all), 'XTickLabel', namesLocComb(iLocComb_all));
ylabel('corr(derived, ideal)');
title('Template reconstruction quality');
grid on;

end
