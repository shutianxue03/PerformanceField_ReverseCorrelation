function data_mirror = getMirror(data)

data_mirror = data;

nDims = ndims(data);
nOri = size(data, nDims -1);

indMirror_L = 1:floor(nOri/2);                 % left side
indMirror_R = flip(ceil(nOri/2+1):nOri);  % right side
indMirror_C = ceil(nOri/2);                     % center

switch nDims
    case 2
        mirror_ave = (data(indMirror_L, :) + data(indMirror_R, :))/2; % average of left and right
        data_mirror(indMirror_L, :) = mirror_ave;
        data_mirror(indMirror_R, :) = mirror_ave;
        data_mirror(indMirror_C, :) = data(indMirror_C, :); % center
    case 3
        mirror_ave = (data(:,indMirror_L, :) + data(:,indMirror_R, :))/2; % average of left and right
        data_mirror(:,indMirror_L, :) = mirror_ave;
        data_mirror(:,indMirror_R, :) = mirror_ave;
        data_mirror(:,indMirror_C, :) = data(:,indMirror_C, :); % center
    case 4
        mirror_ave = (data(:,:,indMirror_L, :) + data(:,:,indMirror_R, :))/2; % average of left and right
        data_mirror(:,:,indMirror_L, :) = mirror_ave;
        data_mirror(:,:,indMirror_R, :) = mirror_ave;
        data_mirror(:,:,indMirror_C, :) = data(:,:,indMirror_C, :); % center
    case 5
        mirror_ave = (data(:,:,:,indMirror_L, :) + data(:,:,:,indMirror_R, :))/2; % average of left and right
        data_mirror(:,:,:,indMirror_L, :) = mirror_ave;
        data_mirror(:,:,:,indMirror_R, :) = mirror_ave;
        data_mirror(:,:,:,indMirror_C, :) = data(:,:,:,indMirror_C, :); % center
end

