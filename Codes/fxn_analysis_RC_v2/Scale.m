function output = Scale(input)
% output = Scale(input)
% Perform an affine scaling to put data in range [0-1].

% [minval,maxval] = MatMinMax(input);
minval = min(input(:));
maxval = max(input(:));
output = (input - minval) ./ (maxval-minval);