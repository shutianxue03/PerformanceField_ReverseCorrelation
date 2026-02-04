%% SX_analysis5_NOM_Trialwise.m
% Author: Shutian Xue
% Purpose: audit existence of trial-wise NOM outputs across all conditions
%
% This version scans all expected file paths and prints a missing-file report.

% =============================
%  Scan all expected outputs
% ==============================

nUnfinshed = 0;

for flag_step=flag_step_all
    for iJob = 1:nJob
        fprintf('--- Job %d/%d ---\n', iJob, nJob)

        for iModelA = iModelA_all
            fprintf('  Model A%d\n', iModelA)

            for iLocSingle = iLocSingle_all

                for isubj = 1:nsubj
                    subjName = subjList{isubj};

                    % Folder must match generating scripts
                    if flag_subjIsHuman
                        nameFolder_NOM_save = sprintf('%s/%s/L%d', nameFolder_Data_NOM_Trialwise, subjName, iLocSingle);
                    else
                        nameFolder_NOM_save = sprintf('%s/%s', nameFolder_Data_NOM_Trialwise, subjName);
                    end

                    switch flag_step
                        case 1 % ---- compIV file ----
                            % expected final output
                            nameFile_compIV = sprintf('%s/n%d_J%d_A%d_compIV.mat', nameFolder_NOM_save, nIter, iJob, iModelA);

                            if ~exist(nameFile_compIV,'file')

                                % 1) If the final file is missing, look for progress-report variants
                                %    e.g., n200_J5_A2_compIV_108_59min.mat
                                nameFile_progress = sprintf('n%d_J%d_A%d_compIV_*min.mat', nIter, iJob, iModelA);
                                nameDir_progress = dir(fullfile(nameFolder_NOM_save, nameFile_progress));

                                if ~isempty(nameDir_progress)

                                    % Print a copy/paste command (terminal) for re-running this combo
                                    fprintf('sbatch shell_NOM_compIV.sh %d %d %d %d %d %d\n', isubj, iLocSingle, iModelA, nIter, nJob, iJob);

                                    % Delete the progress reports of the unfinished job
                                    for iReport = 1:numel(nameDir_progress)
                                        delete(fullfile(nameDir_progress(iReport).folder, nameDir_progress(iReport).name));
                                    end

                                    % Keep track of unfinished jobs
                                    nUnfinished = nUnfinished + 1;  % make sure you initialized this before the loops
                                else
                                    fprintf('%s\nALERT: Neither the finished job nor the progress report (of un unfinished job) exists!\n', nameFile_progress)
                                end
                            end

                        case 2 % ---- fitNOM files (depend on ModelB) ----
                            for iModelB = iModelB_all

                                % expected final output
                                nameFile_fitNOM = sprintf('%s/n%d_J%d_A%dB%d.mat', nameFolder_NOM_save, nIter, iJob, iModelA, iModelB);


                                if ~exist(nameFile_fitNOM,'file')

                                    % 1) If the final file is missing, look for progress-report variants
                                    nameFile_progress = sprintf('n%d_J%d_A%dB%d_compIV_*min.mat', nIter, iJob, iModelA, iModelB);
                                    nameDir_progress = dir(fullfile(nameFolder_NOM_save, nameFile_progress));

                                    if ~isempty(nameDir_progress)

                                        % Print a copy/paste command (terminal) for re-running this combo
                                        fprintf('sbatch shell_NOM_fitNOM.sh %d %d %d %d %d %d %d\n', isubj, iLocSingle, iModelA, iModelB, nIter, nJob, iJob);

                                        % Delete the progress reports of the unfinished job
                                        for iReport = 1:numel(nameDir_progress)
                                            delete(fullfile(nameDir_progress(iReport).folder, nameDir_progress(iReport).name));
                                        end

                                        % Keep track of unfinished jobs
                                        nUnfinished = nUnfinished + 1;  % make sure you initialized this before the loops
                                    else
                                        fprintf('%sL%d %s\nALERT: Neither the finished job nor the progress report (of un unfinished job) exists!\n', subjName, iLocSingle, nameFile_progress)
                                    end
                                end
                            end % iModelB
                    end
                end % isubj
            end % iLocSingle
        end % iModelA
    end % iJob


    fprintf('\n\nNumber of unfinished files: %d\n\n', nUnfinshed)
end