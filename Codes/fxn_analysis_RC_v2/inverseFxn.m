% function x=inverseFxn(y, ifamily, params, x_all)
% dev = nan(length(x_all), 1);
% for ii = 1:length(x_all)
%     x_p = x_all(ii);
% %     y_sim = fxn(x_p, params);
%     y_sim = predSFkernel(x_p, ifamily, params, 0);
%     dev(ii) = abs(y-y_sim);
% end
% [~, ix] = min(dev);
% x = x_all(ix);