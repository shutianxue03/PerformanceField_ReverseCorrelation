function S = GetGrandSeed(nIterPerJob, iJob, nJob, nameFolder)

% -------------------------
% Input checks
% -------------------------
if ~(isnumeric(nIterPerJob) && isscalar(nIterPerJob) && nIterPerJob >= 1)
    error('nIterPerJob must be a positive scalar.');
end
if ~(isnumeric(iJob) && isscalar(iJob) && iJob >= 1)
    error('iJob must be a positive scalar.');
end
if ~(isnumeric(nJob) && isscalar(nJob) && nJob >= 1)
    error('nJob must be a positive scalar.');
end
if nargin < 4 || isempty(nameFolder)
    error('nameFolder is required (directory where GrandSeed.txt will be saved).');
end
nameFolder = char(string(nameFolder));

% Ensure folder exists
if ~exist(nameFolder, 'dir')
    mkdir(nameFolder);
end

% -------------------------
% Create/load grand seed once (race-safe-ish)
% -------------------------
nameFile_GrandSeed = sprintf('%s/GrandSeed.txt', nameFolder);
grandSeed = NaN;

if exist(nameFile_GrandSeed, 'file')
    txt = strtrim(string(fileread(nameFile_GrandSeed)));
    grandSeed = str2double(txt);
else
    grandSeed_candidate = randi(2^31-2);

    fid = -1;
    try
        fid = fopen(nameFile_GrandSeed, 'x'); % atomic create if supported
    catch
        fid = -1;
    end

    if fid >= 0
        fprintf(fid, '%d\n', grandSeed_candidate);
        fclose(fid);
        grandSeed = grandSeed_candidate;
    else
        % Fallback: if file still doesn't exist, create it (not atomic, but prevents deadlock)
        if ~exist(nameFile_GrandSeed, 'file')
            fid2 = fopen(nameFile_GrandSeed, 'w');
            fprintf(fid2, '%d\n', grandSeed_candidate);
            fclose(fid2);
        end

        % Read with retries (handles "file exists but empty" window)
        for k = 1:50
            txt = strtrim(string(fileread(nameFile_GrandSeed)));
            grandSeed = str2double(txt);
            if isfinite(grandSeed), break; end
            pause(0.05);
        end
    end
end

if ~isfinite(grandSeed)
    error('Cannot parse grand seed from: %s', nameFile_GrandSeed);
end

% -------------------------
% Global iteration indices for this job
% -------------------------
iterStart   = (iJob-1)*nIterPerJob + 1;
iterEnd     = min(iJob*nIterPerJob, nJob*nIterPerJob);
iterIdxList = iterStart:iterEnd;

% -------------------------
% Output
% -------------------------
S = struct();
S.grandSeed = grandSeed;
S.iterStart = iterStart;
S.iterEnd   = iterEnd;
S.iterIdxList = iterIdxList;
S.nameFile_GrandSeed = nameFile_GrandSeed;

end
