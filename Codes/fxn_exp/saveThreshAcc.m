
folderName = ['Data/', subjName, '/'];
threshFileName = sprintf('%s%s_ThreshAccRecord.mat', folderName, subjName);
loc_list = params.design.loc_list;

if exp_mode ~= 2 % not practice
    switch exp_mode
        case 0 % titration
            extractThresh
            prev0 = dir(sprintf('%s%s_stair*', folderName, subjName));
%             if isempty(prev0)
                threshT_allSess = nan(1, 5); 
                accT_allSess = nan(1, 5); 
                accExp_allSess = nan(1, 5);
%             end
            threshT_allSess(1,:) = thresh_exp';
            
        case 3 % qstair
            extractThresh
            load(threshFileName)
            threshT_allSess = [threshT_allSess; thresh_exp'];
            if sum(find(loc_list == 2)), accT_allSess = [accT_allSess; record.accPerBlock(1:3), nan, record.accPerBlock(end)];
            else, accT_allSess = [accT_allSess; record.accPerBlock(1), nan, record.accPerBlock(2:end)];
            end
            
        case 1 % exp
            load(threshFileName)
            if isnan(accExp_allSess(1,1))
                if sum(find(loc_list == 2)), accExp_allSess(1,:) = [record.accPerBlock(1:3), nan, record.accPerBlock(end)];
                else, accExp_allSess(1,:) = [record.accPerBlock(1), nan, record.accPerBlock(2:end)];
                end
            else
                if sum(find(loc_list == 2)), accExp_allSess  = [accExp_allSess; [record.accPerBlock(1:3), nan, record.accPerBlock(end)]];
                else, accExp_allSess  = [accExp_allSess; [record.accPerBlock(1), nan, record.accPerBlock(2:end)]];
                end
            end
    end
    save(threshFileName, 'threshT_allSess', 'accT_allSess', 'accExp_allSess')
end