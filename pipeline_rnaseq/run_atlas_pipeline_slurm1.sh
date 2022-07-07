#!/bin/bash
#SBATCH --nodes=1                       # Number of nodes
#SBATCH --ntasks-per-node=1       # number of tasks per node. Equivalent to ppn in TORQUE
#SBATCH --time=30-00:00:00                # Wall time
#SBATCH --partition=p_sslow              # partition. In slurm it is mandatory to set the partition/queue
#SBATCH --mem=3G                      # requested memory per node ( see also : --mem-per-cpu )
#SBATCH --job-name=snakemake          # name of your job. 

module load Anaconda3
source activate atlas_base_kalk2020

snakemake -s atlas_pipeline_slurm1.smk \
          -j 25 \
		  --use-conda \
          --configfile atlas_config.yaml \
          --cluster-config atlas_config_kluster_slurm.yaml \
		  --cluster-status ./status_slurm3.py \
          --cluster "sbatch -t {cluster.time} --job-name={cluster.job-name} --tasks-per-node={cluster.ntasks} --cpus-per-task={cluster.cpus-per-task} --mem-per-cpu={cluster.mem-per-cpu} -p {cluster.partition} -D {cluster.working-directory} --parsable" $1 $2