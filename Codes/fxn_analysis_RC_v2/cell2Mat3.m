

function dataMat3 = cell2Mat3(dataCell)

nlayers = length(dataCell);
[sz1,sz2] = size(dataCell{1});

dataMat3 = nan(sz1, sz2, nlayers);
for n=1:nlayers
    dataMat3(:,:,n) = dataCell{n};
end