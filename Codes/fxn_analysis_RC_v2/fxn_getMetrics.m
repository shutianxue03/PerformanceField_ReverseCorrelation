function metrics = fxn_getMetrics(resp, iPRS, respC, contrast, RT)
pYES = mean(resp == 1);
pHit = mean(iPRS == 1 & resp == 1) / mean(iPRS == 1);
pFA = mean(iPRS == 0 & resp == 1) / mean(iPRS == 0);
pC = mean( (iPRS == 0 & resp == 0) | (iPRS == 1 & resp == 1) );

[dprime, criterion] = SX_sim06_SDT(pHit, pFA);

metrics = [dprime, criterion, [pC, pHit, pFA], nanmean(respC), pYES, mean(1./contrast), median(RT)];
end