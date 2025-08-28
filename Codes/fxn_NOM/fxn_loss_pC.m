function loss = fxn_loss_pC(criterion_potential, pC_titrate, IV_noisy_sim_allT, iPRS_allT)
resp_allT = IV_noisy_sim_allT > criterion_potential;
pHit = sum((iPRS_allT == 1) & (resp_allT == 1)) / sum(iPRS_allT == 1);
pFA = sum((iPRS_allT == 0) & (resp_allT == 1)) / sum(iPRS_allT == 0);
pC_fit = (pHit + 1 - pFA) / 2;
loss = (pC_fit-pC_titrate)^2;
end