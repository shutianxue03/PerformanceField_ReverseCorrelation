function Z = SX_RC_basisProject(e3D_allT, axis_ori_deg, axis_sf_log2, basisOpts)
% SX_RC_basisProject
% Project ORI x SF energy maps into the tensor-product basis used by SX_RC_selectBasis_cv.
% Basis creation is delegated to predSFkernel utility mode.

Bori = predSFkernel('make_basis_ori', axis_ori_deg, basisOpts, 0);
Bsf = predSFkernel('make_basis_sf', axis_sf_log2, basisOpts, 0);

nTrials = size(e3D_allT, 1);
Kori = size(Bori, 2);
Ksf = size(Bsf, 2);
Z = zeros(nTrials, Kori * Ksf);

for iTrial = 1:nTrials
    Ei = squeeze(e3D_allT(iTrial, :, :));
    Zi = Bori' * Ei * Bsf;
    Z(iTrial, :) = Zi(:)';
end
end
