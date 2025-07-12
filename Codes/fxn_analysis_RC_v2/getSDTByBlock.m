

nHit = sum(resp_B & iPrs_B);
nFA = sum(resp_B & (1-iPrs_B));

[dprime_B, criterion_B] = SX_sim06_SDT(nHit, nFA, nTrialsPerBlock/2);

dprime_perSess_perLoc(ib_byLoc(iLoc_B), iLoc_B) = dprime_B;
criterion_perSess_perLoc(ib_byLoc(iLoc_B), iLoc_B) = criterion_B;