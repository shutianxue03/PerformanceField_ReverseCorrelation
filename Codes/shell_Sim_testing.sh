#!/bin/bash -le

#SBATCH --job-name=OODsim
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --output=Logs/OODsim_%A_%a.out
#SBATCH --error=Logs/OODsim_%A_%a.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=vivanxuest@gmail.com
#SBATCH --account=torch_pr_503_general
#SBATCH --array=0-1075

set -euo pipefail

module purge
module load matlab/2025b

cd "/scratch/sx712/PF_RC/Codes"

start=$((SLURM_ARRARY_TASK_ID*5+1))
end=$((start+5))
if [[ ${end} -gt 5376 ]]; then end=5376; fi

for((i=${start}; i<${end}; i++)); do
    line=$(head -${i} /scratch/sx712/PF_RC/Data_job_params/OOD_sim_params.tsv | head -${i} | tail -1)
    IFS=$'\t' read -r \
       noiseCST \
       gaborCST \
       nTrials \
       Nmul_true \
       Nadd_true \
       Nshared_true \
       cSDT_true \
       lambda_whiten \
       iModelB_sim \
       nIter <<< "${line}"
    matlab -batch "try, OOD_sim(${noiseCST}, ${gaborCST}, ${nTrials}, ${Nmul_true}, ${Nadd_true}, ${Nshared_true}, ${cSDT_true}, ${iModelB_sim}, ${nIter}); catch ME, disp(getReport(ME,'extended')); exit(1); end; exit(0);"
done
