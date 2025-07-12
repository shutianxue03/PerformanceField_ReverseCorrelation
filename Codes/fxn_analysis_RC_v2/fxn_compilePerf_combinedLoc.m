
ipA = 1; % obtain pA from 1=all trials; 2=PRS; 3=ABS
ipC = 1;

iLocComb_all_all = {[6,5,3], [6,7], [5,3], [1,8]};
% iLocComb_all_all = {[6,5,3], [6,7], [5,3]};
fprintf('\n\n Compiling behav data: \n')

for ii = iLocComb_all_all
    
    iLocComb_all_ = ii{1};
    
    if length(iLocComb_all_)==3, nameFileLoc = 'all3'; else, nameFileLoc = sprintf('%d%d', iLocComb_all_); end
    fprintf('    L%s...', nameFileLoc)
    
    nameFilePerf_selected = sprintf('%s/n%d_B%d_perf_L%s.mat', nameFolderCompile_BEHAV, nsubj, nB, nameFileLoc);
    dirPerf_selected = dir(nameFilePerf_selected);
    
    for iperf = 1:6 % 1=CS, 2=pA, 3=dprime, 4=criterion, 5=RT, 6=pC
        perf_ave_allSubj = nan(nsubj, nB, length(iLocComb_all_));
        for isubj = 1:nsubj
            subjName = subjList{isubj};
            nSess = nblocks_allSubj(isubj)/5;
            load(sprintf('Data_OOD/%s%d/%s_behavMeas', subjName, nblocks_allSubj(isubj), subjName), '*perSess_perLoc')
            assert(nSess == size(cst_perSess_perLoc, 1));
            if nSess_s<0, indSess = nSess + nSess_s+1:nSess; % take the last few blocks
            else, indSess = nSess_s:nSess; % discard the first few blocks
            end
            if isnan(nSess_s), indSess = 1:nSess; end
            switch iperf
                case 1, perf_perLoc = 1./cst_perSess_perLoc(indSess, :);
                case 2, perf_perLoc = squeeze(pA3_perSess_perLoc(indSess,:, ipA));
                case 3, perf_perLoc = squeeze(dprime_perSess_perLoc(indSess, :));
                case 4, perf_perLoc = squeeze(criterion_perSess_perLoc(indSess, :));
                case 5, perf_perLoc = squeeze(RT_perSess_perLoc(indSess, :));
                case 6, perf_perLoc = squeeze(pC3_perSess_perLoc(indSess, :, ipC));
            end
            
            switch length(iLocComb_all_)
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                case 2
                    switch iLocComb_all_(1)
                        case 1 % 1 vs. 8
                            perf_ave1 = mean(perf_perLoc(:, 1));
                            perf_ave2 = perf_perLoc(:, 2:5); perf_ave2 = mean(perf_ave2(:));
                        case 5 % 5 vs. 3
                            perf_ave1 = mean(perf_perLoc(:, 5));
                            perf_ave2 = mean(perf_perLoc(:, 3));
                        case 6 % 6 vs. 7
                            perf_ave1 = perf_perLoc(:, [2,4]); perf_ave1 = mean(perf_ave1(:));
                            perf_ave2 = perf_perLoc(:, [5,3]); perf_ave2 = mean(perf_ave2(:));
                        case 2 % 2 vs. 4
                            perf_ave1 = mean(perf_perLoc(:, 2));
                            perf_ave2 = mean(perf_perLoc(:, 4));
                        case 4 % 4 vs. 7
                            perf_ave1 = mean(perf_perLoc(:, 4));
                            perf_ave2 = perf_perLoc(:, [5,3]); perf_ave2 = mean(perf_ave2(:));
                    end
                    perf_ave_allSubj(isubj, :, :) = repmat([perf_ave1, perf_ave2], nB, 1);
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                case 3
                    perf_ave6 = perf_perLoc(:, [2,4]); perf_ave6 = mean(perf_ave6(:));
                    perf_ave5 = mean(perf_perLoc(:, 5));
                    perf_ave3 = mean(perf_perLoc(:, 3));
                    perf_ave_allSubj(isubj, :, :) = repmat([perf_ave6, perf_ave5, perf_ave3], nB, 1);
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                case 4
                    perf_ave2 = mean(perf_perLoc(:, 2));
                    perf_ave3 = mean(perf_perLoc(:, 3));
                    perf_ave4 = mean(perf_perLoc(:, 4));
                    perf_ave5 = mean(perf_perLoc(:, 5));
                    perf_ave_allSubj(isubj, :, :) = repmat([perf_ave2, perf_ave4, perf_ave5, perf_ave3], nB, 1);
            end
            
        end % isubj
        
        %%%%%%%%%%%%%
        switch iperf
            case 1, cs_allSubj = perf_ave_allSubj;
            case 2, pA_allSubj = perf_ave_allSubj;
            case 3, dprime_allSubj = perf_ave_allSubj;
            case 4, criterion_allSubj = perf_ave_allSubj;
            case 5, RT_allSubj = perf_ave_allSubj;
            case 6, pC_allSubj = perf_ave_allSubj;
        end
    end % iperf
    
    save(nameFilePerf_selected, 'cs_allSubj', 'pA_allSubj', 'dprime_allSubj', 'criterion_allSubj', 'RT_allSubj', 'pC_allSubj')
    clear *perSess_perLoc
    fprintf(' SAVED\n')
end



