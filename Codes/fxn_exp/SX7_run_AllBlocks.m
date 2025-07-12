function record = SX7_run_AllBlocks(exp_mode, fileName, record, EL, params, demoFlag)

% originally written by Antoine Barbot, adapted by Shutian Xue
% replaced all 'scr.main' to 'screen.wPtr'
% global constant scr visual participant sequence stimulus response timing params

%%
design = params.design;
nTrialsPerBlock = design.nTrialsPerBlock;
% record.iblock_current = record.iblock_start;

for idx_block = 1: design.nBlocks
    % create stim in advance for exp mode
    if exp_mode == 1, exp_createStim4AllTrials, record.run_allTrials = run_allTrials; end
    
    iblock_current = record.iblock_current;
    WaitSecs(.1);
    record.iblock_current = record.iblock_current + 1;
    
    time_blockStart = GetSecs;
    
    SX8_startBlock(params, EL, iblock_current)
   
    run.itrial  = 1; % index of the trial to be saved
%     lastRotInd = 0;
    
    while run.itrial <= nTrialsPerBlock
        %  fprintf('Block %d/%d, trial %d/%d.\n',iblock, design.nBlocks, run.itrial, design.nTrialsPerBlock)
        
        run = SX9_initRun(run.itrial, idx_block, design, record, exp_mode);
        run.ON = EL.ON;
        if ~run.ON, run.fixation = 1; end % if not using eyetracker switch back to 1 in order to save data
        
        %% run one trial
        run = SX10_run_OneTrial(exp_mode, EL, run, params, demoFlag); 
%         if isfield(run, 'lastRotInd'), lastRotInd = run.lastRotInd; else, run.lastRotInd = 0;end
        
        %% update
        if run.fixation && run.check
            % update staircase
            if exp_mode == 1, run.stair = [];
            else, run.stair = usePalamedesStaircase(run.stair,run.correctness);
            end
            
            % save data
            exp_saveData

            % update trial number
            if run.ON, Eyelink('Message', 'TRIAL_END'); end
            run.itrial = run.itrial +1;
        
        else % if fixation is broken, continue to the next trial and put the current trial at the end of the block
            indReorder = [1: run.itrial-1,  run.itrial + 1 : nTrialsPerBlock,  run.itrial];
            if exp_mode ==1, record.run_allTrials = record.run_allTrials(indReorder);
            else, params.design.trialInd = params.design.trialInd(indReorder, :);
            end
        end
    end
    
    record.duration(iblock_current) = GetSecs - time_blockStart;
    record = SX11_endBlock(exp_mode, fileName, params, run, record);
end


