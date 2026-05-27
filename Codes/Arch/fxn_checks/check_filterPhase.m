

% check the phase plot (aper_sz x aper_sz) of the filter
% atan2: measure the angle between two vectors

for ifilterOri = 1:5:nfiltersOri
    for ifilterSF = 1:5:nfiltersSF
        figure
        imagesc(atan2(filter_sin{ifilterOri, ifilterSF}, filter_cos{ifilterOri, ifilterSF}))
        title(sprintf('ORI = %d, SF = %.1f', filtersOri_all(ifilterOri), filtersSF_all(ifilterSF)))
    end
end