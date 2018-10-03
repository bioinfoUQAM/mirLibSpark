#!/bin/bash
#SBATCH --account=def-banire
#SBATCH --time=01:00:00
#SBATCH --nodes=4
#SBATCH --mem=4000M
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=multi-node-spark-181002
#SBATCH --output=jobout/%x-%j.out

module load nixpkgs/16.09
module load spark/2.3.0
#module load python/2.7.14
#module load java/1.8.0_121


export SPARK_IDENT_STRING=$SLURM_JOBID
export SPARK_WORKER_DIR=$SLURM_TMPDIR
start-master.sh

(
export SPARK_NO_DAEMONIZE=1;
srun -x $(hostname -s) -n $((SLURM_NTASKS -1)) --label --output=$SPARK_LOG_DIR/spark-$SPARK_IDENT_STRING-workers.out \
             start-slave.sh -m ${SLURM_MEM_PER_NODE}M -c ${SLURM_CPUS_PER_TASK} spark://$(hostname -f):7077
) &

spark-submit --executor-memory ${SLURM_MEM_PER_NODE}M $SPARK_HOME/pi.py 100000

stop-master.sh


#= sbatch pyspark_submit.sh
#= squeue -u cjwu
