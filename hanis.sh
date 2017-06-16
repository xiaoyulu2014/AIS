#!/bin/bash
#note anything after #SBATCH is a command
#SBATCH --mail-user=xiaoyu.lu@stats.ox.ac.uk
#Email you if job starts, completed or failed
#SBATCH --mail-type=ALL
#SBATCH --job-name=anis_1D
#SBATCH --partition=small 
#Choose your partition depending on your requirements
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#Cpus required for each job
#SBATCH --time=00:40:00
#SBATCH --mem-per-cpu=10000
#Memory per cpu in megabytes          
#SBATCH --array=0-199

# now run your script. Choose one of the following:
matlab -nodisplay -r "hanis; exit"

# print environment variables: the job ID, sub-job's task ID, and sub-job’s job ID.
echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

# Make a directory for output (.txt) and results (e.g. .jld, .json ,.mat) if it doesn't already exist
mkdir -p /data/ziz/not-backed-up/xlu/results/anis/anis_1D_${SLURM_ARRAY_JOB_ID}

# Move experiment outputs & results to the directories made above
mv /data/localhost/not-backed-up/xlu/results/anis/anis_1D_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.mat /data/ziz/not-backed-up/xlu/results/anis/anis_1D_${SLURM_ARRAY_JOB_ID}/${SLURM_ARRAY_TASK_ID}.mat