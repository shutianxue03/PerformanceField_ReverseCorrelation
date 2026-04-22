% Script: simPlot_allCond_compIV_fitNOM.m
%
% Purpose:
% Batch-plot simulation outputs for all IO conditions by loading:
%   (1) compIV outputs
%   (2) fitNOM outputs
% and calling:
%   - NOMplot_compIV
%   - NOMplot_fitNOM
%
% This version matches current saving logic in:
%   - OOD_sim
%   - OOD_NOM_compIV
%   - OOD_NOM_fitNOM
%
% Created by Shutian Xue
% Revised on 2026-04-16

clear; clc; close all;
set(0, 'DefaultFigureVisible', 'off');

%--------------%
SX_RC1_setting;
%--------------%

%% ---------------- settings ----------------
str_part = 'FlexBasis';     % change if needed
iModelA_all = 1;        % or 1:2
iModelB_all = 3;        % or 1:4

nameFolder_Data = sprintf('%s/Data_%s', nameFolder_server, str_part);
nameFolder_Data_OOD = sprintf('%s/Data_OOD_%d%d', nameFolder_Data, nORI, nSF);
nameFolder_Data_NOM_Trialwise = sprintf('%s/Data_NOM_Trialwise_%d%d', nameFolder_Data, nORI, nSF);

nameFolder_Figures_part = fullfile(nameFolder_Figures, sprintf('IO_%s', str_part));
if ~exist(nameFolder_Figures_part, 'dir')
    mkdir(nameFolder_Figures_part);
end


namesMetrics = {'pYES', 'pC', 'pA'};
nMetrics = numel(namesMetrics);

%% ---------------- find IO folders ----------------
nameDir = dir(nameFolder_Data_OOD);
nameDir = nameDir([nameDir.isdir]);
nameDir = nameDir(~ismember({nameDir.name}, {'.', '..'}));
nFiles = numel(nameDir);

if nFiles == 0
    error('No IO folders found in %s', nameFolder_Data_OOD);
end

fprintf('\nFound %d IO folders.\n', nFiles);

%% ---------------- main loop ----------------
for iModelA_fit = iModelA_all
    fprintf('\n\n==================== ModelA = %d ====================\n', iModelA_fit);

    str_loadCompIV = sprintf('n*_A%d_compIV.mat', iModelA_fit);

    for iModelB_fit = iModelB_all
        fprintf('\n------------ ModelB fit = %d (%s) ------------\n', ...
            iModelB_fit, namesModelB{iModelB_fit});

        str_loadfitNOM = sprintf('n*_A%dB%d.mat', iModelA_fit, iModelB_fit);

        for iFile = 1:nFiles

            %% ---------- folder / file names ----------
            nameIO = nameDir(iFile).name;

            nameFolder_OOD_load = fullfile(nameFolder_Data_OOD, nameIO);
            nameFolder_NOM_load = fullfile(nameFolder_Data_NOM_Trialwise, nameIO);

            nameFile_truth = fullfile(nameFolder_OOD_load, 'truth.mat');

            if ~exist(nameFile_truth, 'file')
                fprintf('Missing truth.mat: %s\n', nameFile_truth);
                % continue
            end

            nameDir_compIV = dir(fullfile(nameFolder_NOM_load, str_loadCompIV));
            nameDir_compIV = nameDir_compIV(~contains({nameDir_compIV.name}, 'min'));
            if isempty(nameDir_compIV)
                fprintf('No compIV file found in %s\n', nameFolder_NOM_load);
                % continue
            end
            nameFile_compIV = fullfile(nameDir_compIV(1).folder, nameDir_compIV(1).name);

            nameDir_fitNOM = dir(fullfile(nameFolder_NOM_load, str_loadfitNOM));
            nameDir_fitNOM = nameDir_fitNOM(~contains({nameDir_fitNOM.name}, 'min'));
            if isempty(nameDir_fitNOM)
                fprintf('No fitNOM file found in %s\n', nameFolder_NOM_load);
                % continue
            end
            nameFile_fitNOM = fullfile(nameDir_fitNOM(1).folder, nameDir_fitNOM(1).name);

            %% ---------- progress ----------
            fprintf('\n%s: %d/%d  %s  |  A=%d  B=%d (%s)\n', ...
                datetime('now'), iFile, nFiles, nameIO, ...
                iModelA_fit, iModelB_fit, namesModelB{iModelB_fit});

            %% ---------- parse condition info ----------
            info = fxn_parse_nameIO(nameIO);

            %% ---------- load into structs ----------
            S_truth  = load(nameFile_truth);
            S_compIV = load(nameFile_compIV);
            S_fitNOM = load(nameFile_fitNOM);

            %% ---------- assign loaded fields into workspace ----------
            assign_struct_fields_to_caller(S_truth);
            assign_struct_fields_to_caller(S_compIV);
            assign_struct_fields_to_caller(S_fitNOM);

            %% ---------- variables expected by plotting scripts ----------
            iLocComb = 1;
            subjName = nameIO;
            iModelA = iModelA_fit;

            % Optional condition info for titles / debug
            if ~isempty(fieldnames(info))
                info_this = info; %#ok<NASGU>
            end

            % figure folder: separate by condition + A/B so nothing overwrites
            nameFolder_Figures_perSubj = fullfile(nameFolder_Figures_part, subjName, sprintf('A%dB%d', iModelA_fit, iModelB_fit));
            if ~exist(nameFolder_Figures_perSubj, 'dir')
                mkdir(nameFolder_Figures_perSubj);
            end

            %% ---------- derive convenience variables only if missing ----------
            % Some plotting scripts may expect these.
            if exist('pred_metrics_allIter', 'var') && ~isempty(pred_metrics_allIter)
                pred_test = pred_metrics_allIter{1}; %#ok<NASGU>
            end

            if exist('namesModelBparams', 'var')
                nParams = numel(namesModelBparams{iModelB_fit}); %#ok<NASGU>
            elseif exist('namesModelBparams_short', 'var')
                nParams = numel(namesModelBparams_short{iModelB_fit}); %#ok<NASGU>
            end

            % Only create these if the required variables exist
            if exist('data_train_allIter', 'var') && ~isempty(data_train_allIter)
                if isfield(data_train_allIter{1}, 'IV')
                    NOMc_lb = min(data_train_allIter{1}.IV); %#ok<NASGU>
                    NOMc_ub = max(data_train_allIter{1}.IV); %#ok<NASGU>
                    NOMc_0  = mean([NOMc_lb, NOMc_ub]); %#ok<NASGU>
                end
            end

            %% ---------- Part 1: compIV ----------
            fprintf('%s: Plotting compIV\n', datetime('now'));
            % try
            NOMplot_compIV
            % catch ME
            %     fprintf(2, 'Error in NOMplot_compIV for %s, A=%d, B=%d\n', ...
            %         nameIO, iModelA_fit, iModelB_fit);
            %     fprintf(2, '%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
            % end

            %% ---------- Part 2: fitNOM ----------
            fprintf('%s: Plotting fitNOM\n', datetime('now'));
            % Use the min and max of DV to constrain criterion
            NOMc_lb = min(data_train_allIter{1}.IV);
            NOMc_ub = max(data_train_allIter{1}.IV);
            NOMc_0 = mean([NOMc_lb, NOMc_ub]);
            switch iModelB_fit
                case 1  % FullModel: multi. + additive + shared, criterion
                    params0   = [NOMp1_0,   NOMp2_0,   NOMp3_0, NOMc_0];
                    params_lb = [NOMp1_lb, NOMp2_lb, NOMp3_lb, NOMc_lb];
                    params_ub = [NOMp1_ub, NOMp2_ub, NOMp3_ub, NOMc_ub];

                case 2  % No Nshared: multi. + additive (no shared noise)
                    params0   = [NOMp1_0,   NOMp2_0, NOMc_0];
                    params_lb = [NOMp1_lb, NOMp2_lb, NOMc_lb];
                    params_ub = [NOMp1_ub, NOMp2_ub, NOMc_ub];

                case 3  % No Nmulti: additive + shared (no multi. noise)
                    params0   = [NOMp2_0,   NOMp3_0, NOMc_0];
                    params_lb = [NOMp2_lb, NOMp3_lb, NOMc_lb];
                    params_ub = [NOMp2_ub, NOMp3_ub, NOMc_ub];

                case 4  % No Nadd: multi. + shared (no additive noise)
                    params0   = [NOMp1_0,   NOMp3_0, NOMc_0];
                    params_lb = [NOMp1_lb, NOMp3_lb, NOMc_lb];
                    params_ub = [NOMp1_ub, NOMp3_ub, NOMc_ub];

                case 5  % Just Nadd: additive only
                    params0   = [NOMp2_0, NOMc_0];
                    params_lb = [NOMp2_lb, NOMc_lb];
                    params_ub = [NOMp2_ub, NOMc_ub];

                case 6  % Just Nmul: multi. only
                    params0   = [NOMp1_0, NOMc_0];
                    params_lb = [NOMp1_lb, NOMc_lb];
                    params_ub = [NOMp1_ub, NOMc_ub];

                case 7  % Just Nshared: shared only
                    params0   = [NOMp3_0, NOMc_0];
                    params_lb = [NOMp3_lb, NOMc_lb];
                    params_ub = [NOMp3_ub, NOMc_ub];

                otherwise
                    error('OOD_NOM_Trialwise_fitNOM: Unknown iModelB_fit = %d', iModelB_fit);
            end

            % try
            NOMplot_fitNOM
            % catch ME
            %     fprintf(2, 'Error in NOMplot_fitNOM for %s, A=%d, B=%d\n', ...
            %         nameIO, iModelA_fit, iModelB_fit);
            %     fprintf(2, '%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
            % end

        end % iFile
    end % iModelB_fit
end % iModelA_fit

fprintf('\nDone.\n');

%% Local helper functions
function assign_struct_fields_to_caller(S)
fn = fieldnames(S);
for i = 1:numel(fn)
    assignin('caller', fn{i}, S.(fn{i}));
end
end

function info = fxn_parse_nameIO(nameIO)

info = struct();

tok = regexp(nameIO, ...
    'cG([-\d\.eE]+).*Nm([-\d\.eE]+).*Na([-\d\.eE]+).*Ns([-\d\.eE]+).*Cz([-\d\.eE]+).*B(\d+)', ...
    'tokens', 'once');

assert(~isempty(tok), 'ALERT: Could not parse nameIO: %s', nameIO);

info.gaborCST     = str2double(tok{1})/100;
info.Nmul_true    = str2double(tok{2});
info.Nadd_true    = str2double(tok{3});
info.Nshared_true = str2double(tok{4});
info.Cz_true      = str2double(tok{5});
info.iModelB_sim  = str2double(tok{6});
end

