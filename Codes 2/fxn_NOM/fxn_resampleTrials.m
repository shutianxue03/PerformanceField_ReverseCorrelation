
e3D_train = [];
iPRS_train = [];
resp_train = [];
cst_train = [];
respC_train = [];

e3D_test = [];
iPRS_test = [];
resp_test = [];
cst_test = [];
respC_test = [];
iPair_test=[];

for iiLoc=1:length(iLoc_all)

    dataMtx_nonrand = [];
    e3D_nonrand_train = [];
    e3D_nonrand_test = [];
    iPRS_nonrand = [];
    iPair_nonrand = [];
    resp_nonrand = [];
    cst_nonrand = [];

    for iSS = 1:length(iSess_select)
        indLoc = (dataMatrix(:, 5) == iLoc_all(iiLoc)) & (dataMatrix(:, 2) == iSess_select(iSS));
        dataMtx_nonrand = [dataMtx_nonrand; dataMatrix(indLoc, :)];
        e3D_nonrand_train = cat(1, e3D_nonrand_train, e3D_allT_train(indLoc, :, :));
        e3D_nonrand_test = cat(1, e3D_nonrand_test, e3D_allT_test(indLoc, :, :));
        iPRS_nonrand = [iPRS_nonrand; dataMatrix(indLoc, 6)];
        iPair_nonrand = [iPair_nonrand; dataMatrix(indLoc, 8)];
        resp_nonrand = [resp_nonrand; dataMatrix(indLoc, 9)];
        cst_nonrand = [cst_nonrand; dataMatrix(indLoc, 11)];
    end % iSS
    
    cst_unik_nonrand = unique(cst_nonrand);
    iPair_unik_nonrand = unique(iPair_nonrand); % range: [1, ntrials_ofThisSubj/2]
    
    %% sample trials data from each SINGLE location according to the unik pair
    rng('shuffle')
    indUnikPair_rand = randperm(ntrials_perSingleLoc/2, ntrials_perSingleLoc/2); % range: [1, ntrials_perSingleLoc/2]
    
    nPairs_train = ntrain_perSingle/2;
    nPairs_test = ntrials_perSingleLoc/2 - nPairs_train; assert(nPairs_train+nPairs_test == length(indUnikPair_rand))
    indUnikPair_rand_train_ = indUnikPair_rand(1:nPairs_train); % range: [1, ntrialsPerLoc/2]
    indUnikPair_rand_test_ = indUnikPair_rand(nPairs_train+1:end); assert(length(indUnikPair_rand_test_) == nPairs_test) % range: [1, ntrialsPerLoc/2]
    
    % TRAIN group
    [e3D_train_, iPRS_train_, resp_train_, cst_train_, respC_train_] = fxn_getRespC(e3D_nonrand_train, indUnikPair_rand_train_, ...
        iPair_nonrand, iPair_unik_nonrand, iPRS_nonrand, resp_nonrand, cst_nonrand);
    
    % TEST group
    [e3D_test_, iPRS_test_, resp_test_, cst_test_, respC_test_, iPair_test_] = fxn_getRespC(e3D_nonrand_test, indUnikPair_rand_test_, ...
        iPair_nonrand, iPair_unik_nonrand, iPRS_nonrand, resp_nonrand, cst_nonrand);
    
    e3D_train = cat(1, e3D_train, e3D_train_);
    iPRS_train = [iPRS_train; iPRS_train_];
    resp_train = [resp_train; resp_train_];
    cst_train = [cst_train; cst_train_];
    respC_train = [respC_train; respC_train_];
    
    e3D_test = cat(1, e3D_test, e3D_test_);
    iPRS_test = [iPRS_test; iPRS_test_];
    resp_test = [resp_test; resp_test_];
    cst_test = [cst_test; cst_test_];
    iPair_test = [iPair_test; iPair_test_];
    respC_test = [respC_test; respC_test_];
    
end % iLoc_
