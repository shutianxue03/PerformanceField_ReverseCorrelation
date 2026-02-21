%% ========================= helper: partial eta^2 for 2-way RM =========================
function [eta2p_A, eta2p_B, eta2p_AB] = fxn_eta2p_rm2(Y)
% Y: [nSubj x nA x nB] (here: A=ModelB, B=Loc)
% Returns partial eta^2 for A, B, and AxB.
%
% Assumes reasonably complete data; uses 'omitnan' for safety.

[nSubj, nA, nB] = size(Y);

% Means
grand = mean(Y, [1 2 3], 'omitnan');              % scalar
S  = mean(Y, [2 3], 'omitnan');                   % [nSubj x 1 x 1]
A  = mean(Y, [1 3], 'omitnan');                   % [1 x nA x 1]
B  = mean(Y, [1 2], 'omitnan');                   % [1 x 1 x nB]
AB = mean(Y, 1, 'omitnan');                       % [1 x nA x nB]
SA = mean(Y, 3, 'omitnan');                       % [nSubj x nA x 1]
SB = mean(Y, 2, 'omitnan');                       % [nSubj x 1 x nB]

% SS for effects
SS_A  = nB*nSubj * sum((squeeze(A) - grand).^2, 'omitnan');
SS_B  = nA*nSubj * sum((squeeze(B) - grand).^2, 'omitnan');

tmpAB = AB - A - B + grand;                        % [1 x nA x nB]
SS_AB = nSubj * sum(tmpAB(:).^2, 'omitnan');

% SS for corresponding error terms
tmpSA = SA - S - A + grand;                        % [nSubj x nA x 1]
SS_SA = nB * sum(tmpSA(:).^2, 'omitnan');

tmpSB = SB - S - B + grand;                        % [nSubj x 1 x nB]
SS_SB = nA * sum(tmpSB(:).^2, 'omitnan');

% Residual for interaction (subject x A x B)
tmpSAB = Y - SA - SB - AB + A + B + S - grand;     % [nSubj x nA x nB]
SS_SAB = sum(tmpSAB(:).^2, 'omitnan');

% Partial eta^2
eta2p_A  = SS_A  / (SS_A  + SS_SA);
eta2p_B  = SS_B  / (SS_B  + SS_SB);
eta2p_AB = SS_AB / (SS_AB + SS_SAB);
end