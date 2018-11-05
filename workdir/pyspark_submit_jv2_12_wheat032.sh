#!/bin/bash
#SBATCH --job-name=jv2-12-wheat-mirL11500M1n32pn-181028
#SBATCH --account=def-banire
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --mem=115000M
#SBATCH --cpus-per-task=32
#SBATCH --ntasks-per-node=1
#SBATCH --error=jobout/%x-%j.err
#SBATCH --output=jobout/%x-%j.out
#SBATCH --mail-user=g39103001@gm.ym.edu.tw
#SBATCH --mail-type=ALL

#= maximun --cpus-per-task=32
#= maximun --mem=115000M

#= preloaded: python2.7, perl, java
#= module loaded: pyspark, duskmasker, bowtie, RNAfold
#= included dependencies: miranda, VARNA
module load nixpkgs/16.09
module load spark/2.3.0
module load gcc/5.4.0
module load viennarna/2.4.9
module load bowtie/1.1.2
module load blast+/2.6.0

#= python requirements: statsmodels (with this, it includes: numpy, scipy, pandas, patsy), seaborn (with this, it includes: matplotlib)
pip install --user requests
pip install --user -r requirements.txt

export _JAVA_OPTIONS="-Xms8g -Xmx16g"
export SPARK_IDENT_STRING=$SLURM_JOBID
export SPARK_WORKER_DIR=$SLURM_TMPDIR
start-all.sh

sleep 5
MASTER_URL=$(grep -Po '(?=spark://).*' $SPARK_LOG_DIR/spark-${SPARK_IDENT_STRING}-org.apache.spark.deploy.master*.out)

NWORKERS=$((SLURM_NTASKS - 1))
SPARK_NO_DAEMONIZE=1 
srun -n ${NWORKERS} -N ${NWORKERS} --label --output=$SPARK_LOG_DIR/spark-%j-workers.out start-slave.sh -m ${SLURM_MEM_PER_NODE}M -c ${SLURM_CPUS_PER_TASK} ${MASTER_URL} &
slaves_pid=$!


#= example:
#spark-submit --master ${MASTER_URL} --executor-memory ${SLURM_MEM_PER_NODE}M /home/cjwu/project/cjwu/gitRepo/mirLibSpark/workdir/pi.py 1000

spark-submit --master ${MASTER_URL} --executor-memory ${SLURM_MEM_PER_NODE}M ../src/mirLibPipeline.py --species wheat --input_type fasta --adapter TGGAATTCTCGGGTGCCAAGGAACTC


kill $slaves_pid
stop-all.sh



#= sbatch pyspark_submit_jv2_02.sh
#= squeue -u cjwu
#= scancel <jobid>

