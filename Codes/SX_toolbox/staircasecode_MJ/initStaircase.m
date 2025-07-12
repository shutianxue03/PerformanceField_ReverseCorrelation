function staircase = initStaircase(stairParams)
% Initialize the staircase
switch stairParams.whichStair
    case 1 % best PEST
        % uniform prior with the mode selected as xCurrent
        meanMode = 'mode';
        prior = ones(1,length(stairParams.alphaRange));
        prior = prior/sum(prior); % make a uniform prior
    case 2 % QUEST
        % normally-distributed prior with the mean selected as xCurrent
        meanMode = 'mean';
        prior = PAL_pdfNormal(stairParams.alphaRange,stairParams.questMean,stairParams.questSD);
end

% set prior to be last posterior if it is provided
if isfield(stairParams,'lastPosterior') && ~isempty(stairParams.lastPosterior)
    prior = stairParams.lastPosterior;
end

staircase = PAL_AMRF_setupRF('priorAlphaRange',stairParams.alphaRange,...
    'stopcriterion','trials','stoprule',inf,'beta',stairParams.fitBeta,...
    'lambda',stairParams.fitLambda,'gamma',stairParams.fitGamma,...
    'meanmode',meanMode,'PF',stairParams.PF,'prior',prior);
staircase.threshPerformance = stairParams.threshPerformance;
staircase.type = 'maxlikelihood';
staircase.updateAfterTrial = stairParams.updateAfterTrial;
staircase.preUpdateLevels = stairParams.preUpdateLevels;
if staircase.updateAfterTrial > 0, staircase.xCurrent = staircase.preUpdateLevels(1); end
