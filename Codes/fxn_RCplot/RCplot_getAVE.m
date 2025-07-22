
%% empty containers
marg_medB_allSubj_all = cell(2, ntypes, nLoc8);
marg_CI_B_neg_allSubj_all = marg_medB_allSubj_all;
marg_CI_B_pos_allSubj_all = marg_medB_allSubj_all;
marg_aveSubj_all = marg_medB_allSubj_all;
marg_semSubj_all = marg_medB_allSubj_all;

pred_medB_allSubj_all = cell(2, ntypes,  ncomb6);
pred_CI_B_neg_allSubj_all = pred_medB_allSubj_all;
pred_CI_B_pos_allSubj_all = pred_medB_allSubj_all;
pred_aveSubj_all = pred_medB_allSubj_all;
pred_semSubj_all = pred_medB_allSubj_all;

params_medB_allSubj_all = cell(2, ntypes, max([nparamsOri, nparamsSF]), ncomb6);
params_CI_B_neg_allSubj_all = params_medB_allSubj_all;
params_CI_B_pos_allSubj_all = params_medB_allSubj_all;
params_aveSubj_all = params_medB_allSubj_all;
params_semSubj_all = params_medB_allSubj_all;

R2_medB_allSubj_all = cell(2, ntypes,  ncomb6);
R2_CI_B_neg_allSubj_all = R2_medB_allSubj_all;
R2_CI_B_pos_allSubj_all = R2_medB_allSubj_all;
R2_aveSubj_all = R2_medB_allSubj_all;
R2_semSubj_all = R2_medB_allSubj_all;


%%
for ifeature = 1:2
    switch ifeature
        case 1
            nfilters = nfiltersOri;
            nparams_full = nparamsOri;
            marg_allSubj = margORI_allSubj; % nsubj x nB x ntypes x nfilters x nLoc
            margPred_allSubj = margPred_ORI_allSubj; % nsubj x nB x ntypes x ncomb x nlines x nfilters
            margParams_allSubj = margParams_ORI_allSubj; % nsubj x nB x ntypes x ncomb x nlines x nparams_full
            margR2_allSubj= margR2_ORI_allSubj; % nsubj x nB x ntypes x ncomb x nlines x 2 (slope and intercept)
        case 2
            nfilters = nfiltersSF;
            nparams_full = nparamsSF;
            marg_allSubj = margSF_allSubj;
            margPred_allSubj = margPred_SF_allSubj;
            margParams_allSubj = margParams_SF_allSubj;
            margR2_allSubj= margR2_SF_allSubj;
    end
    
    % reorder the location (to select the dominant eye side for HM)
%     marg_allSubj = reorderLoc(marg_allSubj, eyeDom);
%     margPred_allSubj = reorderLoc(margPred_allSubj, eyeDom);
%     margParams_allSubj = reorderLoc(margParams_allSubj, eyeDom);
%     margR2_allSubj = reorderLoc(margR2_allSubj, eyeDom);

    for itype = 1:3
        for icomb = 1:ncomb6
            %%%%%%
            %   marg  %
            %%%%%%
            for iline = 1:2
                iLoc = combInd(icomb, iline);
               
                [medB_allSubj, CI_B_neg_allSubj, CI_B_pos_allSubj, aveSubj, semSubj] = getIDVD_GROUP_metrics(marg_allSubj(:, :,  itype, :, iLoc), ntrialsProp);
                marg_medB_allSubj_all{ifeature, itype, iLoc} = medB_allSubj;
                marg_CI_B_neg_allSubj_all{ifeature, itype, iLoc} = CI_B_neg_allSubj;
                marg_CI_B_pos_allSubj_all{ifeature, itype, iLoc} = CI_B_pos_allSubj;
                marg_aveSubj_all{ifeature, itype, iLoc} = aveSubj;
                marg_semSubj_all{ifeature, itype, iLoc} = semSubj;
            end
            
            %%%%%%%%%
            %   margPred  %
            %%%%%%%%%
            [medB_allSubj, CI_B_neg_allSubj, CI_B_pos_allSubj, aveSubj, semSubj] = getIDVD_GROUP_metrics(margPred_allSubj(:, :,  itype, icomb, :, :), ntrialsProp);
            pred_medB_allSubj_all{ifeature, itype, icomb} = medB_allSubj;
            pred_CI_B_neg_allSubj_all{ifeature, itype, icomb} = CI_B_neg_allSubj;
            pred_CI_B_pos_allSubj_all{ifeature, itype, icomb} = CI_B_pos_allSubj;
            pred_aveSubj_all{ifeature, itype, icomb} = aveSubj;
            pred_semSubj_all{ifeature, itype, icomb} = semSubj;
            
            %%%%%%%%
            %   params  %
            %%%%%%%%
            for iparam = 1:nparams_full
                [medB_allSubj, CI_B_neg_allSubj, CI_B_pos_allSubj, aveSubj, semSubj] = getIDVD_GROUP_metrics(margParams_allSubj(:, :,  itype, icomb, :, iparam), nan);
                params_medB_allSubj_all{ifeature, itype, iparam, icomb} = medB_allSubj;
                params_CI_B_neg_allSubj_all{ifeature, itype, iparam, icomb} = CI_B_neg_allSubj;
                params_CI_B_pos_allSubj_all{ifeature, itype, iparam, icomb} = CI_B_pos_allSubj;
                params_aveSubj_all{ifeature, itype, iparam, icomb} = aveSubj;
                params_semSubj_all{ifeature, itype, iparam, icomb} = semSubj;
            end % end of iparam
            
            %%%%%%%%
            %   marg R2  %
            %%%%%%%%
            [medB_allSubj, CI_B_neg_allSubj, CI_B_pos_allSubj, aveSubj, semSubj] = getIDVD_GROUP_metrics(margR2_allSubj(:, :,  itype, icomb, :), nan);
            R2_medB_allSubj_all{ifeature, itype, icomb} = medB_allSubj;
            R2_CI_B_neg_allSubj_all{ifeature, itype, icomb} = CI_B_neg_allSubj;
            R2_CI_B_pos_allSubj_all{ifeature, itype, icomb} = CI_B_pos_allSubj;
            R2_aveSubj_all{ifeature, itype, icomb} = aveSubj;
            R2_semSubj_all{ifeature, itype, icomb} = semSubj;
            
        end % end of icomb
    end % end of itype
end % end of ifeature

%% helper fxns
%%%%%%%%%%%%%%%%%%%
function [medB_allSubj, CI_B_neg_allSubj, CI_B_pos_allSubj, aveSubj, semSubj] = getIDVD_GROUP_metrics(data, ntrialsProp)

% get median and CI 
[medB_allSubj, ~, ~, CI_B_neg_allSubj, CI_B_pos_allSubj] = getCI(data, 1, 2, 0);

% get group AVE and SEM
if isnan(ntrialsProp)
    [aveSubj, ~, ~, semSubj, ~] = getCI(medB_allSubj, 2, 1, 0);
else
    [aveSubj, ~, ~, semSubj, ~] = getCI(medB_allSubj, 2, 1, 0, ntrialsProp);
end

% squeeze the median
medB_allSubj = squeeze(medB_allSubj);
CI_B_neg_allSubj = squeeze(CI_B_neg_allSubj);
CI_B_pos_allSubj = squeeze(CI_B_pos_allSubj);
end


