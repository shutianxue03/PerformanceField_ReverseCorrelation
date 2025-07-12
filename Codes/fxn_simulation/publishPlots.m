for n = 10:5:35
    load(['fxn_simulation/dataMat_nCST',num2str(n)])
    plot_sim
end
makeBeep(ones(1,3) * .2, [1e3, 750, 500])