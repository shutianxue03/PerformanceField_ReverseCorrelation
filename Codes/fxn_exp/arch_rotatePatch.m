% function [patch_r, rand1]  = rotatePatch(patch, s_, lastRotInd)
% rand1_ = randperm(8) - 1;
% rand1 = rand1_(1);

% while rand1 == lastRotInd % make sure the toration is not the same in two consecutive trials
%     rand1_ = randperm(8) - 1;
%     rand1 = rand1_(1);
% end

% patch_r = imrotate(patch, rand1 * 45); % rotated by 0, 90, 180 or 270 deg

% % cut the edge
% s_r = size(patch_r , 1);
% start_ = round((s_r- s_)/2) + 1;
% ind = start_ : (start_ + s_ -1 );
% if length(ind) ~= s_
%     error('rotated by %d deg, the size of the cut image (%d) does not match the Gabor (%d)\n', rand1(1) * 45, length(d), s_)
% end
% patch_r = patch_r(ind,ind);