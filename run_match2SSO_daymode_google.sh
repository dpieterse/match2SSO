#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --mem=0
#SBATCH --job-name=m2sso-day

#SBATCH --partition=p1gb4t

#SBATCH --open-mode=append
#SBATCH --output=/home/sa_105685508700717199458/RunMatch2SSO/log/Slurm/daymode_%A.log
#SBATCH --mail-user=danielle@blackgem.org
#SBATCH --mail-type=ALL
##SBATCH --mail-type=FAIL,TIME_LIMIT

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Ask the OS for a free ephemeral port
CLUSTER_PROXY_PORT=$(python3 -c '
import socket
s = socket.socket()
s.bind(("", 0))
print(s.getsockname()[1])
s.close()
')
export CLUSTER_PROXY_PORT

# Start up a tunnel to the Slurm login node, so that the downloading of MPCORB
# will use the whitelisted IP address of that node.
ssh -D "$CLUSTER_PROXY_PORT" -N -f "$SLURM_SUBMIT_HOST"

# Run match2SSO
/opt/apps/singularity/3.11.0/bin/singularity exec /home/sa_105685508700717199458/Containers/MLBG_latest.sif python /Software/match2SSO/match2SSO.py --mode day --telescope BG4

# Close the ssh connection
pkill -f "ssh -D $CLUSTER_PROXY_PORT"
