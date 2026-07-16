#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=4GB
#SBATCH --job-name=m2sso-day

#SBATCH --account=b19-meerlicht-ag
#SBATCH --reservation=meerlicht
#SBATCH --partition=Main
#SBATCH --output=/idia/projects/meerlicht/RunMatch2SSO/log/Slurm/daymode_%A.log
#SBATCH --mail-user=d.pieterse@astro.ru.nl
#SBATCH --mail-type=FAIL,TIME_LIMIT,REQUEUE

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

/software/common/singularity/4.4.1/bin/singularity exec /idia/projects/meerlicht/Containers/MLBG_latest.sif python /Software/match2SSO/match2SSO.py --mode day
