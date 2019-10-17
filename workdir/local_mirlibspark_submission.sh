#!/bin/bash
#SBATCH --job-name=localspark.module
#SBATCH --time=00:10:00
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=110G
#
#SBATCH --account=your-account
#SBATCH --mail-user=your@email.com
#SBATCH --mail-type=ALL
#SBATCH --error=%x-%j.err
#SBATCH --output=%x-%j.out

## --------------------------------------
## author: Chao-Jung Wu
## --------------------------------------


#= preloaded: python2.7, perl, java
#= module loaded: pyspark, duskmasker, bowtie, RNAfold
#= included dependencies: miranda, VARNA
#= viennarna will load python3.7, so need to reload python2.7.14 to override
module load nixpkgs/16.09
module load gcc/7.3.0
module load viennarna/2.4.11
module load bowtie/1.1.2
module load blast+/2.7.1
module load spark/2.4.4
module load python/2.7.14

#= python requirements: statsmodels (includes: numpy, scipy, pandas, patsy), 
#                       seaborn (includes: matplotlib)
virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate
pip install --no-index --upgrade pip
pip install --no-index -r requirements.txt


#= JAVA memory allocation (space for else than RDD operations)
export _JAVA_OPTIONS="-Xms3g -Xmx10g"
#export SPARK_IDENT_STRING=$SLURM_JOBID
#export SPARK_WORKER_DIR=$SLURM_TMPDIR


spark-submit --master local[*] \
             --executor-memory ${SLURM_MEM_PER_NODE}M \
             ../src/mirLibPipeline.py \
             --jobid ${SLURM_JOBID}
