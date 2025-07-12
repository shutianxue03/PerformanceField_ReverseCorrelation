function [fileName,record] = SX4_checkFile(exp_mode, subjName, params)

design = params.design;
stairParams = params.stairParams;

%% fileName
folderName = sprintf('Data/%s', subjName);
switch exp_mode
    case 0, exp_name = 'stair';
    case 1, exp_name = 'exp';
    case 2, exp_name = 'practice';
    case 3, exp_name = 'qstair';
end
fileName.brief = sprintf('%s/%s_%s', folderName, subjName, exp_name);
fileName.brief_dir = dir([fileName.brief,'*']); % find out all possible files

%% generate data container
blank = nan(design.nBlocks * 30, design.nTrialsPerBlock);

%% load titration data
if ~isempty(fileName.brief_dir), fileExist = 1; else, fileExist = 0; end
switch exp_mode
    case 0 % titration mode
        % keep doing titration
        if fileExist 
            if length(fileName.brief_dir) > 1 % just in case not all loc were finished last time
                fileName.formal = fileName.brief_dir(end).name;
            else
                fileName.formal = fileName.brief_dir.name;
            end
            load(fileName.formal)
            nblocks = sum(~isnan(record.correctness (:,1)));
            if ~isfloat(nblocks / design.nCuedLoc), extractThresh, end % to show the cst thresh at all loci
            
            % adjust prior based on the pdf from data collected so far
            design = record.params.design;
            load(sprintf('Data/%s/%s_PDF.mat', subjName, subjName), 'lastPosterior_allLoc')
            for iLoc = 1:design.nCuedLoc
                for istair = 1:design.nStairs
                    record.stairs_all{iLoc,istair}.prior = lastPosterior_allLoc(iLoc, :);
                end
            end
        % just start titration
        else 
            mkdir(folderName) % create a new folder for the subject
            stairs_all = SX5_initStaircases(design, stairParams, exp_mode);
            record = exp_createRecord(params, blank, stairs_all);
        end
        
    case 1 % exp
        thresh_exp = input('\n\n         >>> Enter the threshold to use: ')
        
        if fileExist
            load(fileName.brief_dir(end).name)
        else
            stairs_all = SX5_initStaircases(design, stairParams, exp_mode);
            record = exp_createRecord(params, blank, stairs_all);
            record.iblock_current = 1;
        end
        record.thresh_exp = thresh_exp;
        
    case 2 % practice
        stairs_all = SX5_initStaircases(design, stairParams, exp_mode);
        record = exp_createRecord(params, blank, stairs_all);
        
    case 3 % qstair
        % load the saved PDF file: lastPosterior_allLoc, x, mean, and std
        load([folderName, subjName, '_PDF'], 'lastPosterior_allLoc', 'x')
        stairParams.lastPosterior_allLoc = lastPosterior_allLoc;
        stairParams.alphaRange = x;
        
        % define new staircases
        stairs_all = SX5_initStaircases(design, stairParams, exp_mode);
        record = exp_createRecord(params, blank, stairs_all);
        if fileExist
            num = regexp(fileName.brief_dir(end).name, 'B[\d]+', 'match');
            num = convertStringsToChars(num{1});
            record.iblock_current = str2double(num(2))+1;
        else, record.iblock_current = 1;
        end
end

fprintf('Current progress: finished %d blocks.\n', record.iblock_current-1)

