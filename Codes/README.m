// %% README

// % NOTE: running the experimental & analysis scripts requires installation
// % of psychtoolbox (the version with which scripts were developed is 3.0.18)

// % [Local]: codes are run on the local computer; parfor loops might be
// % involved
// % [OOD]: computationally expensive processes will be run on HPC (high-performance computer) maintained by
// % NYU; thus scripts and data transfer are required (refer to xx for a list of orders)

// %% Abbreviations:
// %     ORI: orientation
// %     SF: spatial frequency
// %     PRS: present (i.e., Gabor-present)
// %     ABS: absent (i.e., Gabor-absent)
// %     cst: contrast
// %     pC: proportion of correctness, or accuracy
// %     pA: proportion of agreement, or response consistency
// % .   RT: response time 

// %% Navigate the scripts
// %==================================================================
// % Run Model Simulation
// %==================================================================
// % 1. Run OOD_sim_v2.m on HPC to simulate trial-wise stimuli and data, which are saved in Data_OOD_2929 and Data_NOM_Trialwise_2929 on the server
// % If run on HPC, remeber to modify "nameFolder_server" in SX_RC1_setting.m;

// % 2. Download the simulated data from the server to the local machine
// % 3. Analyze simulated data "IO_xx" using SX_analysis5_NOM_trialWise.m


// %==================================================================
// % RUN EXP
// %==================================================================
// % PF_RC_runExp
// % functions cited are in the folder fxn_exp

// % Outputs are saved in /Data/subjName/ 
// % For the titration process: 
// % 1. performances are saved in the file with the name formatted as
// % subjName_stair_B#L#_yearmonthdateTtime(24 hour) overides the staircases
// % saved before 
// % 2. the PDF for running the PEST, file name formatted as subjName_PDF.mat

// % For the main experiment: 
// % 1. patches presented on all trials & performances of EACH block
// % with file name formatted as subjName_exp_B#L#_yearmonthdateTtime

// % e.g., Data/AS/AS_exp_B001L5_20220218T1207.mat
// % corresponds to the first block completed by a subject whose initial is
// % AS, in the main experiment on Feb 18th, 2022 at 12:07, in which
// % the target is located at the Lower Vertical Meridian (indicated by the
// % number after 'L', 1=Fovea, 2=Left, 3=Upper,4=Right, 5=Lower)

// % 2. eye data are saved in /eyedata/ (not transfered from L1 computer yet)

// %==================================================================
// % analyze the titration data (to determine the RMS cst of the noise patch)
// %==================================================================
// % PF_RC_analysis_titration
// % helper functions are saved in /fxn_analysis_titration

// %==================================================================
// % Reverse-correlation pipeline
// %==================================================================
// % helper functions are saved in /fxn_analysis_RC_v2, fxn_MC and fxn_model

// %%%%%%%%%%%%%%%%%
// % Step 1: [local] run SX_analysis1_extractData.m 
// %%%%%%%%%%%%%%%%%
// % to extract behav data and derive source energy (of all trials)
// % Step 1.1. check corr between energy bins
// % fxn_analysis_RC/check_ebins_corr

// % loadAndDelete: for each subj, load a certain file, process saved variables, and save all

// %%%%%%%%%%%%%%%%%
// % Step 2: [Local] run SX_analysis2_kernels.m to do RC analysis on ALL trials (no resampling)
// % code is similar to OOD_boot (noresampling)
// % Purpose: use the sensitivity kernels derived from ALL trials to do model comparison 
// %               (i.e., to determine which model best captures ORI and SF tuning functions) 
// %%%%%%%%%%%%%%%%%

// %%%%%%%%%%%%%%%%%
// % Step 3: [OOD] run shell_all_MC.sh, which submits jobs running OOD_MC
// % [local] SX_analysis3_MC
// % conduct model comparison to decide which fxn to fit ORI/SF kernels

// %%%%%%%%%%%%%%%%%
// % Step 4: [OOD] run shell_all_boot.sh, which submits jobs running OOD_boot
// % [local] run SX_analysis4_Boot to compile the boostrapped kernels/fittings
// % of two locations (all compiled are saved in Data_compile/)
// % [local] comp_params_vs_tuningC: compare the estimated parameters and tuning characteristics

// %%%%%%%%%%%%%%%%%
// % Step 5. Plot and present
// % fxn_analysis_RC_v2/plotAll_Loc2 to produce figures and run stats to compare 2 locations
// % fxn_analysis_RC_v2/plot_IDVD: plot idvd behav meas, assess whether need to exclude any observer
// % VSS2023/code/plotAll_Loc4 to produce figures and run stats to compare 4 or 5 locations (fovea + 4 polar angles)
// % talk2022/code/pre3_plot: make plots for presentation purposes: schematic logistic regression; schematic tuning fxns; schematic gain & tuning
// % codes for plotting corr between NOM params and CS is in SX_analysis5_model
