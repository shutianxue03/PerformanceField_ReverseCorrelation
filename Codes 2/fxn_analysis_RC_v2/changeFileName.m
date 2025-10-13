% changeFileName
subjList = {'SX', 'SP', 'LS', 'RE'};
% subjList = {'LS'};
nsubj = length(subjList);

for isubj = 1:nsubj
    subjName = subjList{isubj};
    folderName = ['Data/', subjName, '/'];
    dataFileDir = dir([folderName, subjName, '_exp*']); % find out all possible files
    
    nblocks = length(dataFileDir);
    for iblock = 1:nblocks
        nameOld = dataFileDir(iblock).name;
        if nameOld(9) == '0' % 01-09
            nameNew = [nameOld(1:8), '0', nameOld(9:end)];
        elseif (nameOld(9) == '1') && (nameOld(11) == 'L') % 10
            nameNew = [nameOld(1:8), '0', nameOld(9:end)];
        elseif (nameOld(9) == '1') && (nameOld(10) ~= '0') % 11-19
            nameNew = [nameOld(1:9), '0', nameOld(10:end)];
        elseif nameOld(9) > 1 % 20-99
            nameNew = [nameOld(1:8), '0', nameOld(9:end)];
        end

        movefile([folderName, nameOld], [folderName, nameNew])
    end
end