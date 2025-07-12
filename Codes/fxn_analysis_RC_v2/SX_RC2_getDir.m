
nameFolder =  sprintf('Data_OOD/%s%d/', subjName, nblocks);
if isempty(dir(nameFolder)), mkdir(nameFolder), end

nameBehavMeas = sprintf('Data_OOD/%s%d/%s_behavMeas.mat', subjName, nblocks, subjName);
dirBehavMeas = dir(nameBehavMeas);

nameTgtPatch = sprintf('Data_OOD/%s%d/%s_patch.mat', subjName, nblocks, subjName);
dirTgtPatch = dir(nameTgtPatch);

nameNoisePatch = sprintf('Data/%s/%s%d_noiseP.mat', subjName, subjName, nblocks);
dirNoisePatch = dir(nameNoisePatch);

namePatchOutputs = sprintf('Data/%s/%s_patch_outputs.mat', subjName, subjName);
dirPatchOutputs = dir(namePatchOutputs);

%%
nameEnergy = sprintf('Data_OOD/%s%d_energy_%d_%d.mat', subjName, nblocks, nfiltersOri, nfiltersSF);
dirEnergy = dir(nameEnergy);

% nameIDVDdata = sprintf('Data_OOD/%s%d_B%d_%d_%d.mat', subjName, nblocks, nB, nfiltersOri, nfiltersSF);
% dirIDVDdata = dir(nameIDVDdata);

%%
namePublishIDVD = sprintf('publishedPDFs/%s%d_B%d_RC.pdf', subjName, nblocks, nB);
dirPublishIDVD = dir(namePublishIDVD);

