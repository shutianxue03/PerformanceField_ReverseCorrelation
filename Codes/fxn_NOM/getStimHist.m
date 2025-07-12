
% given that the BG contrast has been changed, only look at trials that
% take bigger proportion
iLoc = 1;
pp = patch_both_perComb_p{iLoc};
np = length(pp);

energy_PRS = nan(np/2, 1);
energy_ABS = energy_PRS;
for ip = 1:(np/2)
    pp_PRS = pp{ip};
    pp_ABS = pp{ip+np/2};
    %     figure, imshow(pp{ip})
    if pp_PRS(1,1) > .4, energy_PRS(ip) = sumsqr(pp_PRS)/std(pp_PRS(:)); end
    if pp_ABS(1,1) > .4, energy_ABS(ip) = sumsqr(pp_ABS)/std(pp_ABS(:)); end
    
    
end
    
%     sumsqr(pp{ip})
%     std(pp{ip}(:))
%     energy_PRS(ip) = sumsqr(pp{ip})/std(pp{ip}(:));
%     energy_ABS_ = sumsqr(pp{ip+np/2});
%     if energy_PRS_>2000
% %         ind_high = [ind_high, ip];
%         energy_PRS_high = [energy_PRS_high, energy_PRS_];
%         energy_ABS_high = [energy_ABS_high, energy_ABS_];
%     else
%         energy_PRS_low = [energy_PRS_low, energy_PRS_];
% %         ind_low = [ind_low, ip];
%         energy_ABS_low= [energy_ABS_low, energy_ABS_];
    
    
% end

%%
figure, hold on
histogram(energy_PRS, 100,'FaceColor', 'r')
% histogram(energy_PRS_low, 100,'FaceColor', 'r')
histogram(energy_ABS, 100, 'FaceColor', 'b')
b')
% stem(ind_low, 'b')