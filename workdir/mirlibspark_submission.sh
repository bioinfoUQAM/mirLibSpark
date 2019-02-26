#!/bin/bash
#SBATCH --job-name=mirLibSpark
#SBATCH --account=youraccount
#SBATCH --mail-user=yourname@email.com
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --mem=20000M
#SBATCH --cpus-per-task=32
#SBATCH --ntasks-per-node=1
#SBATCH --error=%x-%j.err
#SBATCH --output=%x-%j.out
#SBATCH --mail-type=ALL


## --------------------------------------
## 0. Preparation
## --------------------------------------
#= preloaded: python2.7, perl, java
#= module loaded: pyspark, duskmasker, bowtie, RNAfold
#= included dependencies: miranda, VARNA
#= viennarna will load python3.7, so need to reload python2.7.14 to override
module load nixpkgs/16.09
module load spark/2.3.0
module load gcc/5.4.0
module load viennarna/2.4.9
module load bowtie/1.1.2
module load blast+/2.6.0
module load python/2.7.14
#
#= python requirements: statsmodels (includes: numpy, scipy, pandas, patsy), 
#                       seaborn (includes: matplotlib)
pip install --user requests
pip install --user -r requirements.txt
#
#= JAVA memory allocation (space for else than RDD operations)
export _JAVA_OPTIONS="-Xms3g -Xmx6g"
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
NWORKERS=$((SLURM_NTASKS - 1))
export SPARK_NO_DAEMONIZE=1 
srun -n ${NWORKERS} -N ${NWORKERS} --label --output=$SPARK_LOG_DIR/spark-%j-workers.out start-slave.sh -m ${SLURM_MEM_PER_NODE}M -c ${SLURM_CPUS_PER_TASK} ${MASTER_URL} &
slaves_pid=$!



## --------------------------------------
## 3. Submit a task to the Spark cluster
## --------------------------------------
spark-submit --master ${MASTER_URL} --executor-memory ${SLURM_MEM_PER_NODE}M ../src/mirLibPipeline.py

## --------------------------------------
## 4. Clean up
## --------------------------------------
kill $slaves_pid
stop-all.sh



## --------------------------------------
## footnotes
## --------------------------------------
#= https://www.sherlock.stanford.edu/docs/software/using/spark/

