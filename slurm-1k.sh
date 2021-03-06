#!/bin/sh
##SBATCH --account=pi-josewalt
#SBATCH  --clusters=faculty
#SBATCH  --partition=isecc
#SBATCH  --qos=isecc
#SBATCH --time=24:00:00
#SBATCH  --nodes=1
#SBATCH  --ntasks-per-node=12
#SBATCH --mem=120000
#SBATCH  --array=1-1000
#SBATCH  --job-name="1kCSP"
#SBATCH  --output="./dat/console/1kCSP-%a.out"
#SBATCH  --error="./dat/console/1kCSP-%a.err"
##SBATCH --mail-user=caigao@buffalo.edu
##SBATCH --mail-type=ALL
#SBATCH --exclude=cpn-p26-[15-18]
##SBATCH --requeue

##echo "SLURM_JOB_ID="$SLURM_JOB_ID
##echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
##echo "SLURM_NNODES"=$SLURM_NNODES
##echo "SLURMTMPDIR="$SLURMTMPDIR
##echo "SLURM_ARRAYID="$SLURM_ARRAYID
##echo "SLURM_ARRAY_JOB_ID"=$SLURM_ARRAY_JOB_ID
##echo "SLURM_ARRAY_TASK_ID"=$SLURM_ARRAY_TASK_ID
##echo "working directory = "$SLURM_SUBMIT_DIR

#NPROCS=`srun --nodes=${SLURM_NNODES}$ bash -c 'hostname' |wc -l`
#echo "NPROCS="$NPROCS

echo "start"

./bin/pulse ${SLURM_ARRAY_TASK_ID}

echo "stop"


