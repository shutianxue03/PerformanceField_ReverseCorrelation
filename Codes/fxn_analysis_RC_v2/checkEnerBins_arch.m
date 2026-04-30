
% function [rHit, rFA, ebinPRS, ebinABS, efficiency] = checkEnerBins(energy, data_both, ifilterORI, ifilterSF, nLoc, nbins_e)

% binsEdge_e = linspace(0, 1, nbins_e+1);
% ebinPRS = nan(nLoc, nbins_e);
% ebinABS = ebinPRS;
% rHit = ebinPRS;
% rFA = ebinPRS;
% efficiency = nan(nLoc, 3); % 3 indicates efficincy, d_real and d_ideal

% %%
% for iLoc = 1:nLoc
%     e = squeeze(energy(iLoc, :, :, :));
    
%     ntrials = size(e,1)/2;
%     e_prs = e(1:ntrials, ifilterORI, ifilterSF);
%     e_abs = e(ntrials+1:end, ifilterORI, ifilterSF);
%     resp_prs = data_both(1:ntrials, iLoc, 2); % the 2nd column is ANSWER (1=YES 0=NO)
%     resp_abs = data_both(ntrials+1:end, iLoc, 2);
    
%     % get efficiency
% %     e_prs_mean = mean(e_prs);
% %     e_abs_mean = mean(e_abs);
% %     e_prs_var = var(e_prs);
% %     e_abs_var = var(e_abs);
% %     d_ideal = (e_prs_mean - e_abs_mean)/sqrt((e_prs_var + e_abs_var)/2);
% %     d_real = norminv(mean(resp_prs))-norminv(mean(resp_abs));
% %     efficiency(iLoc, :) = [(d_real/d_ideal)^2, d_real, d_ideal]; % efficiency
    
%     % get hit and FA rate by energy level
%     edge_prs = quantile(e_prs, binsEdge_e); % min, several medians, and max
%     edge_abs = quantile(e_abs,binsEdge_e);
    
%     ebin_prs_allBins = nan(1, nbins_e-1); 
%     rHit_allBins = ebin_prs_allBins;
%     ebin_abs_allBins = ebin_prs_allBins;      
%     rFA_allBins = ebin_prs_allBins;
    
%     for ibin = 1:nbins_e % nbins_e = length(nbins_e)-1
%         ind_prs = (e_prs >=  edge_prs(ibin)) & (e_prs <=  edge_prs(ibin+1));
%         ebin_prs_allBins(ibin) = mean(e_prs(ind_prs));
%         rHit_allBins(ibin) = mean(resp_prs(ind_prs)); % say YES given gabor PRS
        
%         ind_abs = (e_abs >=  edge_abs(ibin)) & (e_abs <=  edge_abs(ibin+1));
%         ebin_abs_allBins(ibin) = mean(e_abs(ind_abs));
%         rFA_allBins(ibin) = mean(resp_abs(ind_abs)); % say YES given gabor ABS
%     end

%     ebinPRS(iLoc, :)  = ebin_prs_allBins;
%     ebinABS(iLoc, :)  = ebin_abs_allBins;
%     rHit(iLoc, :) = rHit_allBins;
%     rFA(iLoc, :) = rFA_allBins;
% end




