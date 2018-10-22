#!/bin/bash
#SBATCH --job-name=jv2-20-wheat-mirL11500M1n32pn-181022-newsubmitfile
#SBATCH --account=def-banire
#SBATCH --time=00:10:00
#SBATCH --nodes=2
#SBATCH --mem=115000M
#SBATCH --cpus-per-task=16
#SBATCH --ntasks-per-node=2
#SBATCH --error=jobout/%x-%j.err
#SBATCH --output=jobout/%x-%j.out
#SBATCH --mail-user=wu.chaojung@gmail.com
#SBATCH --mail-type=ALL

#= maximun --cpus-per-task=32
#= maximun --mem=115000M



## --------------------------------------
## 0. Preparation
## --------------------------------------
#= preloaded: python2.7, perl, java
#= module loaded: pyspark, duskmasker, bowtie, RNAfold
#= included dependencies: miranda, VARNA
module load nixpkgs/16.09
module load spark/2.3.0
module load gcc/5.4.0
module load viennarna/2.4.9
module load bowtie/1.1.2
module load blast+/2.6.0
#
#= python requirements: statsmodels (with this, it includes: numpy, scipy, pandas, patsy), seaborn (with this, it includes: matplotlib)
pip install --user requests
pip install --user -r requirements.txt
#
#= JAVA memory allocation (space for else than RDD operations)
export _JAVA_OPTIONS="-Xms3g -Xmx4g"
#= identify the Spark cluster with the Slurm jobid
export SPARK_IDENT_STRING=$SLURM_JOBID
export SPARK_WORKER_DIR=$SLURM_TMPDIR



## --------------------------------------
## 1. Start the Spark cluster master
## --------------------------------------
start-all.sh
sleep 5
MASTER_URL=$(grep -Po '(?=spark://).*' $SPARK_LOG_DIR/spark-${SPARK_IDENT_STRING}-org.apache.spark.deploy.master*.out)



## --------------------------------------
## 2. Start the Spark cluster workers
## --------------------------------------
#export SPARK_WORKER_CORES=${SLURM_CPUS_PER_TASK:-1}
#export SPARK_MEM=$(( ${SLURM_MEM_PER_CPU:-4096} * ${SLURM_CPUS_PER_TASK:-1} ))M
#export SPARK_DAEMON_MEMORY=$SPARK_MEM
#export SPARK_WORKER_MEMORY=$SPARK_MEM
#export SPARK_EXECUTOR_MEMORY=$SPARK_MEM
#
## start the workers on each node allocated to the job
#export SPARK_NO_DAEMONIZE=1
#srun  --output=$SPARK_LOG_DIR/spark-%j-workers.out --label \
#      start-slave.sh ${MASTER_URL} &
#
################################################
NWORKERS=$((SLURM_NTASKS - 1))
export SPARK_NO_DAEMONIZE=1 
srun -n ${NWORKERS} -N ${NWORKERS} --label --output=$SPARK_LOG_DIR/spark-%j-workers.out start-slave.sh -m ${SLURM_MEM_PER_NODE}M -c ${SLURM_CPUS_PER_TASK} ${MASTER_URL} &
slaves_pid=$!



## --------------------------------------
## 3. Submit a task to the Spark cluster
## --------------------------------------
#= example:
#spark-submit --master ${MASTER_URL} --executor-memory ${SLURM_MEM_PER_NODE}M /home/cjwu/project/cjwu/gitRepo/mirLibSpark/workdir/pi.py 1000
#
#spark-submit --master ${MASTER_URL} --executor-memory ${SLURM_MEM_PER_NODE}M ../src/mirLibPipeline.py ../paramfile_WHEAT_IWGSC_graham.txt
spark-submit --master ${MASTER_URL} --executor-memory ${SLURM_MEM_PER_NODE}M ../src/mirLibPipeline.py ../paramfile_ATH_TAIR10_graham.txt



## --------------------------------------
## 4. Clean up
## --------------------------------------
kill $slaves_pid
stop-all.sh



## --------------------------------------
## footnotes
## --------------------------------------
#= sbatch pyspark_submit_jv2_02.sh
#= squeue -u cjwu
#= scancel <jobid>
#
#= https://www.sherlock.stanford.edu/docs/software/using/spark/

