#!/bin/bash
#SBATCH --account=def-banire
#SBATCH --time=00:30:00
#SBATCH --nodes=2
#SBATCH --mem=100000M
#SBATCH --cpus-per-task=6
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=jv2-03-mirL10000M2n6pn-181016
#SBATCH --error=jobout/%x-%j.err
#SBATCH --output=jobout/%x-%j.out
#SBATCH --mail-user=wu.chaojung@gmail.com
#SBATCH --mail-type=ALL


#= preloaded: python2.7, perl, java
#= module loaded: pyspark, duskmasker, bowtie, RNAfold
#= included dependencies: miranda, VARNA
module load nixpkgs/16.09
module load spark/2.3.0
module load gcc/5.4.0
module load viennarna/2.4.9
module load bowtie/1.1.2
module load blast+/2.6.0

#= python requirements: statsmodels (includes: numpy), seaborn
pip install --user requests
pip install --user -r requirements.txt

export _JAVA_OPTIONS="-Xms2g -Xmx4g"
export SPARK_IDENT_STRING=$SLURM_JOBID
export SPARK_WORKER_DIR=$SLURM_TMPDIR
start-all.sh

sleep 5
MASTER_URL=$(grep -Po '(?=spark://).*' $SPARK_LOG_DIR/spark-${SPARK_IDENT_STRING}-org.apache.spark.deploy.master*.out)

NWORKERS=$((SLURM_NTASKS - 0))
SPARK_NO_DAEMONIZE=1 srun -n ${NWORKERS} -N ${NWORKERS} --label --output=$SPARK_LOG_DIR/spark-%j-workers.out start-slave.sh -m ${SLURM_MEM_PER_NODE}M -c ${SLURM_CPUS_PER_TASK} ${MASTER_URL} &
slaves_pid=$!


#= example:
#= spark-submit --master ${MASTER_URL} --executor-memory ${SLURM_MEM_PER_NODE}M /home/cjwu/project/cjwu/gitRepo/mirLibSpark/workdir/cedar_training/pi.py 1000
#srun -n ${NWORKERS} -N ${NWORKERS} spark-submit --master ${MASTER_URL} --executor-memory ${SLURM_MEM_PER_NODE}M /home/cjwu/project/cjwu/gitRepo/mirLibSpark/workdir/cedar_training/pi.py 1000
srun -n ${NWORKERS} -N ${NWORKERS} spark-submit --master ${MASTER_URL} --executor-memory ${SLURM_MEM_PER_NODE}M ../src/mirLibPipeline.py ../paramfile_ATH_TAIR10_graham.txt

kill $slaves_pid
stop-all.sh





#= sbatch pyspark_submit_jv2_02.sh
#= squeue -u cjwu
#= scancel <jobid>

