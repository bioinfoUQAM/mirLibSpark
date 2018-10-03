#!/bin/bash
#SBATCH --account=def-banire
#SBATCH --time=01:00:00
#SBATCH --nodes=4
#SBATCH --mem=4000M
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=multi-node-spark-181002
#SBATCH --output=jobout/%x-%j.out
#SBATCH --mail-user=wu.chaojung@gmail.com
#SBATCH --mail-type=ALL

#= preloaded: python2.7, perl, java
#= module loaded: pyspark, duskmasker, bowtie, RNAfold
#= dependencies: miranda, VARNA
module load nixpkgs/16.09
module load spark/2.3.0
module load gcc/5.4.0
module load viennarna/2.4.9
module load bowtie/1.1.2
module load blast+/2.6.0


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

#181002-2031
#Submitted batch job 12281602
#
#Submitted batch job 12281943
