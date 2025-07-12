function seprblty = getSeparability(kernels_raw)
% getSeparability: Computes separability of a 2D kernel matrix
%
% This function assesses the separability of a 2D kernel matrix by
% comparing it with a reconstruction from marginals across spatial
% frequency (SF) and orientation (ORI). The separability score is
% determined by the correlation between the original and reconstructed
% matrices.
%
% INPUT:
%   kernels_raw   - 2D kernel matrix with values across ORI and SF dimensions
%
% OUTPUT:
%   seprblty      - Separability score

% Calculate marginal means across spatial frequency and orientation dimensions
margSF = mean(abs(kernels_raw), 1);  % Mean across orientations (row-wise)
margORI = mean(abs(kernels_raw), 2); % Mean across spatial frequencies (column-wise)

% Reconstruct the kernel matrix assuming separability (outer product of marginals)
kernels_re = margORI * margSF;

% Calculate separability as the correlation between original and reconstructed kernel
seprblty = corr(kernels_raw(:), kernels_re(:));

