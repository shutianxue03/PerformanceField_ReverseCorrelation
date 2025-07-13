function stairParams = validateParams(stairParams)

% below lines are my code to set default, NOT require setting the field
% name

% below is Michael's code, requires setting the field name
% check if all parameters were set, if not, then set to the default
setParams = fieldnames(stairParams);
setParams = sort(setParams);

for i = 1:length(setParams)
    switch setParams{i}
        
        % setting possible threshold estimates
        case 'alphaRange'
            if isempty(stairParams.alphaRange)
                stairParams.alphaRange = 0.01:0.01:1;
                fprintf('ALPHA RANGE: Set to 0.01:0.01:1 (DEFAULT)\n');
            end
            
        % check that inputtted PF is accurate
        case 'PF'
            if ~isfield(stairParams,'PF')
                stairParams.PF = @arbWeibull;
                fprintf('PSYCHOMETRIC function: Set to arbWeibull (DEFAULT)\n');
            else
                if ~isa(stairParams.PF,'function_handle')
                    stairParams.PF = eval(['@',stairParams.PF]);
                elseif ~ismember(func2str(stairParams.PF),{'arbWeibull' 'arbLogistic' 'PAL_Weibull' 'PAL_Quick' 'PAL_Gumbel' ...
                        'PAL_HyperbolicSecant' 'PAL_Logistic'})
                    error('PSYCHOMETRIC function: Inputted PF is not in the possible list of options.');
                end
            end
            
        % check for proper beta
        case 'fitBeta'
            if ~isfield(stairParams,'fitBeta')
                stairParams.fitBeta = 2;
                fprintf('BETA: Set to 2 (DEFAULT)\n');
            end
            
            % check for proper gamma
        case 'fitGamma'
            if ~isfield(stairParams,'fitGamma')
                stairParams.fitGamma = 0.5;
                fprintf('GAMMA: Set to 0.5 (DEFAULT)\n');
            end
            
            % check for proper lambda
        case 'fitLambda'
            if ~isfield(stairParams,'fitLambda')
                stairParams.fitLambda = 0.01;
                fprintf('LAMBDA: Set to 0.01 (DEFAULT)\n');
            end
            
            % check if experimenter wants to continue staircase from previous run
        case 'lastPosterior'
            if isempty(stairParams.lastPosterior)
                stairParams.lastPosterior = [];
            elseif ~isfield(stairParams,'lastPosterior')
                stairParams.lastPosterior = [];
            else
                fprintf('POSTERIOR: Staircase will continue from previously computed posterior.\n');
            end
            
            
            % check if experimenter wants to use an arbritrary threshold performance
        case 'threshPerformance'
            if isempty(stairParams.threshPerformance)
                stairParams.threshPerformance = 0.75;
            end
            % if using arbWeibull, force user to input a threshold
            % performance
            if ~isempty(strfind(func2str(stairParams.PF),'arb')) && isempty(stairParams.threshPerformance)
                error(['EXPECTED PERFORMANCE: If using an arbitrary (modified) performance level you '...
                    'need to input that expected performance level.']);
            end
            
            % choosing staircases
        case 'whichStair'
            if ismember(stairParams.whichStair,[1 2])
                if stairParams.whichStair==2
                    % check for mean and sd input of prior
                    questInput = {'questMean' 'questSD'};
                    questDefault = {'stairParams.alphaRange(round(length(stairParams.alphaRange)/2))' '1'};
                    idx = ismember(questInput,setParams);
                    if any(idx)
                        % check if values are empty
                        questVal = {stairParams.(questInput{1}), stairParams.(questInput{2})};
                        emptyIdx = find(cellfun('isempty',questVal));
                        if ~all(isnumeric([stairParams.questMean stairParams.questSD])) && isempty(emptyIdx)
                            % check that values inputted are numbers
                            error('QUEST: mean and standard deviation should be numbers');
                        else
                            % set empty values to the default
                            for q = 1:length(emptyIdx)
                                fprintf(['Setting ',questInput{emptyIdx(q)}, ' to ', questDefault{emptyIdx(q)}, ' (DEFAULT)\n']);
                                stairParams.(questInput{emptyIdx(q)}) = eval(questDefault{emptyIdx(q)});
                            end
                        end
                    else
                        % set the empty or the non-set input to the default
                        for q = 1:length(questInput), stairParams.(questInput{q}) = eval(questDefault{q}); end
                    end
                end
            else
                fprintf('whichStair: Setting to use best PEST (DEFAULT)\n');
                stairParams.whichStair = 1;
            end
            
            % trial number, after which, staircase will update
        case 'updateAfterTrial'
            if isempty(stairParams.updateAfterTrial) || ~isfield(stairParams,'updateAfterTrial')
                stairParams.updateAfterTrial = 0;
            end
            
            % levels to display before the staircase gets updated
        case 'preUpdateLevels'
            if isempty(stairParams.preUpdateLevels) || ~isfield(stairParams,'preUpdateLevels')
                if stairParams.updateAfterTrial>0
                    stairParams.preUpdateLevels = repmat(median(stairParams.alphaRange),1,stairParams.updateAfterTrial);
                else
                    stairParams.preUpdateLevels = [];
                end
            elseif numel(stairParams.preUpdateLevels) ~= stairParams.updateAfterTrial
                warning('# of pre-update levels need to match number of pre-update trials. Using default (median of alpha range)');
                stairParams.preUpdateLevels = repmat(median(stairParams.alphaRange),1,stairParams.updateAfterTrial);
            end
    end
end
