
function [energy_allBins, pYES_allBins] = checkEnerBins(energy, resp, nbins_e)

binsEdge_e = linspace(0, 1, nbins_e+1);

edge = quantile(energy, binsEdge_e); % min, several medians, and max

energy_allBins = nan(1, nbins_e-1);
pYES_allBins = energy_allBins;

for ibin = 1:nbins_e % nbins_e = length(nbins_e)-1
    ind = (energy >=  edge(ibin)) & (energy <=  edge(ibin+1));
    energy_allBins(ibin) = mean(energy(ind));
    pYES_allBins(ibin) = mean(resp(ind)); % say YES
end


