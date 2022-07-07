#!/bin/bash
#SBATCH --nodes=1                       # Number of nodes
#SBATCH --ntasks-per-node=1       # number of tasks per node. Equivalent to ppn in TORQUE
#SBATCH --time=3-00:00:00                # Wall time
#SBATCH --partition=p_medium              # partition. In slurm it is mandatory to set the partition/queue
#SBATCH --mem=8G                      # requested memory per node ( see also : --mem-per-cpu )
#SBATCH --job-name=snakemake          # name of your job. 

module load Anaconda3
#source activate atlas_base_kalk2020
conda activate /home/debbiarm/.conda/envs/atlas_base_kalk2020/envs/atlas_base_kalk2020_new

snakemake -s atlas_pipeline_slurm2.smk \
          -j 15 \
          --use-conda \
          --conda-frontend conda \
          --configfile atlas_config.yaml \
          --cluster-config atlas_config_kluster_slurm.yaml \
          --cluster-status ./status_slurm.py \
          --max-jobs-per-second 0.1 \
          --max-status-checks-per-second 5 \
          --restart-times 10 \
          --latency-wait 60 \
          --keep-going \
          --scheduler greedy \
          --cluster "sbatch -t {cluster.time} --job-name={cluster.job-name} --tasks-per-node={cluster.ntasks} --cpus-per-task={cluster.cpus-per-task} --mem-per-cpu={cluster.mem-per-cpu} -p {cluster.partition} -D {cluster.working-directory} --parsable" $1 $2 $3