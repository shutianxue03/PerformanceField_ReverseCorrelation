function [dprime, criterion] = SX_sim06_SDT(pHit, pFA)

ntrials = 1e5;

zHit = norminv(pHit);
zFA = norminv(pFA);

if zHit == -Inf, zHit = norminv(1/ntrials); end
if zHit == Inf, zHit = norminv((ntrials-1)/ntrials); end
if zFA == -Inf, zFA = norminv(1/ntrials); end
if zFA == Inf, zFA = norminv((ntrials-1)/ntrials); end

dprime =  zHit - zFA;
criterion = -(zHit + zFA)/2;

