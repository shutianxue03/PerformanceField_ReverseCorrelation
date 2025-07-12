
noiseVersion = input(sprintf('How do you like the noise patch:\n   - 0 = the same noise patch all the time\n   - 1 = same, but rotated\n   - 2 = 2 types of patches\n   - 3 = 3 types of patches'));
noiseVersion = randperm(4,1)-1;
noiseVersion = 1;
% 0 = the same noise patch all the time
% 1 = same, but rotated
% 2 = 2 types of patches
% 3 = 3 types of patches

numNoisePatch = noiseVersion;
if noiseVersion==0, numNoisePatch=1;end

% create and save a noise patch
noisePatch_dir = dir(['Data/noisePatch_',subjName,'.mat']);

if isempty(noisePatch_dir)
    noisePatch = cell(1, numNoisePatch);
    for ip = 1:numNoisePatch

        noiseLooksGood = 0;
        while noiseLooksGood == 0
            noisePatch{ip} = exp_createNoisePatch(noise, stim, ip);
            noiseLooksGood = input(['noise patch #', num2str(ip),' looks great (1) or not (0): ']);
        end
    end
    save(sprintf('Data/noisePatch_%s', subjName), 'noisePatch')
else
    load(['Data/', noisePatch_dir.name], 'noisePatch')
end

noise.noiseVersion = noiseVersion;
noise.numNoisePatch = numNoisePatch;
noise.noisePatch = noisePatch;
fprintf('==============\nTarget SF = %d, , noise cst = %d%%.\n==============\n', gaborSF, round(noise.contrast*100))
