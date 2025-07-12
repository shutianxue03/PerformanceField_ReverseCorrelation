
% 
% %% plot PF: HM, VM, peripheral
% kernels2D_HM = squeeze(mean(kernels2D_ave(:, [2,4] , :, :),2));
% kernels2D_VM = squeeze(mean(kernels2D_ave(:, [3,5] , :, :),2));
% kernels2D_peri = squeeze(mean(kernels2D_ave(:, [2:5] , :, :),2));
% 
% figure('Position', [0 0 length(titles_PF)*300 300])
% for n=1:length(titles_PF)
%     switch n
%         case 1, kernels2D_ = kernels2D_HM;
%         case 2, kernels2D_ = kernels2D_VM;
%         case 3, kernels2D_ = kernels2D_peri;
%     end
%     for itype = 3%1:ntypes
%         % 
%         subplot(1, 3, n), hold on
%         
%         imagesc(filtersOri_all-90, flip(filtersSF_all_log), flip(squeeze(kernels2D_(itype, :, :))'))
%         colorbar%, caxis([min(kernels2D_()), max(max(kernels2D_ave(itype, :,:)))])
%         axis square
%         
%         xlabel('ori')
%         xticks([-90, 0, 90])
%         xticklabels([-90,0,90])
%         ylabel('SF')
%         yticks(filtersSF_all_log([1,round(nfiltersSF/2),nfiltersSF])) % yticks vector must be increasing
%         yticklabels(filtersSF_all([1,round(nfiltersSF/2),nfiltersSF]))
%         title(titles_PF{n})
%     end
% end
% set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)
% sgtitle('PF (average of kernels)')
% 
% 
% %% plot diff between kernels
% % diffType1_all = [ones(1,7), 6, 5, 6, 6]; % center for 6 times, HM, S, HM, HN
% % diffType2_all = [2:8, 7, 3, 3, 5]; % L, R, N, S, HM, VM, peripheral, VM, N, N, S
% diffType1_all = [6,1,5];
% diffType2_all = [7,8,3];
% nComp = length(diffType1_all);
% 
% 
% %%
% figure('Position', [0 200 nComp * 300 300])
% for icomp = 1:nComp
%     diffType1 = diffType1_all(icomp);
%     diffType2 = diffType2_all(icomp);
%     
%     for itype = 3%1:ntypes
%         subplot(ntypes, nComp, icomp + nComp*(itype-1))
%         subplot(1,3,icomp)
%         kernels2D1 = squeeze(kernels2D_ave_comp{diffType1}(itype, :, :));
%         kernels2D2 = squeeze(kernels2D_ave_comp{diffType2}(itype, :, :));
%         
%         imagesc(filtersOri_all-90, flip(filtersSF_all_log), flip((kernels2D1 - kernels2D2)'))
%        
%         colorbar
%         xlabel('ori')
%         xticks([-90, 0, 90])
%         xticklabels([-90,0,90])
%         ylabel('SF')
%         yticks(filtersSF_all_log([1,round(nfiltersSF/2),nfiltersSF])) % yticks vector must be increasing
%         yticklabels(filtersSF_all([1,round(nfiltersSF/2),nfiltersSF]))
%         
%         axis square
%         title(sprintf('%s vs. %s', diffTypeNames{diffType1}, diffTypeNames{diffType2}))
%     end
% end
% set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)
% 
