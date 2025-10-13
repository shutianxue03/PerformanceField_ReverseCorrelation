filtered_noise = exp_CreateFilteredNoise(noise);
limit0to1 = @(x) min(max(x,0),1);

cst_all = [.01, .1, .5, 1, 2];
nCST = length(cst_all);
energy = nan(nCST,nfilters);
figure, hold on
for icst = 1:nCST
    gabor = exp_CreateGabor(stim, cst_all(icst), 0);
    
    stimsPrs_01 = limit0to1(noise.ratio_base+gabor*noise.ratio_gaborInTgt+filtered_noise).*stim.mask + noise.ratio_base*(1-stim.mask);;
    stimsAbs_01 = limit0to1(noise.ratio_base+filtered_noise).*stim.mask + noise.ratio_base*(1-stim.mask);
    
    plot4parts = 0;
    lumiBG = stimsPrs_01(1,1);
    energy_ = nan(1,nfilters);
    for ff = 1:nfilters
        [a,b] = SX_sim04_getEnergy(stimsPrs_01 , lumiBG, SFfilter_sin{ff}, SFfilter_cos{ff}, plot4parts,0, filterSF_all(ff));
        energy_(ff) = a;
    end
    energy(icst, :) = energy_;
    
    plot(energy_)
    legends{icst} = num2str(cst_all(icst));
end


legend(legends)

