function stairs_all = SX5_initStaircases(design, stairParams, exp_mode)

if exp_mode == 0, subjName = input('Enter subjName to set the prior dist (hit ENTER if unnesscary): ');
else, subjName = [];
end

stairs_all = cell(5, design.nStairs); % though fewer than 5 loci are tested, the size of container is 5
for ss = 1:design.nStairs
    for ll = design.loc_list
        % narrow down the alpha range for quick titration
        if exp_mode == 3
            stairParams.lastPosterior = stairParams.lastPosterior_allLoc(ll,:); 
        end

        % premap the alpha space at the beginning of all titration
        if (exp_mode == 0) && (isempty(subjName)) % only premap the alpha space in 'Titration', not 'quick titration'
            if ss == 1, stairParams.preUpdateLevels = stairParams.preUpdateLevels;
            else, stairParams.preUpdateLevels = fliplr(stairParams.preUpdateLevels);
            end
        end
        
        % set up staircase 
        stairs_all{ll,ss} = usePalamedesStaircase(stairParams);
        
        % reconstruct prior given previous data
        if ~isempty(subjName)
            load(sprintf('Data/%s/%s_PDF.mat', subjName, subjName), 'lastPosterior_allLoc')
            stairs_all{ll,ss}.prior = lastPosterior_allLoc(iLoc, :);
        end
    end
end
clear ss ll
