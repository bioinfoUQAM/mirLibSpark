#!/bin/bash
#SBATCH --job-name=mirprefer
#SBATCH --account=def-banire
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --mem=20000M
#SBATCH --cpus-per-task=32
#SBATCH --ntasks-per-node=1
#SBATCH --error=jobout_ath_experiment/%x-%j.err
#SBATCH --output=jobout_ath_experiment/%x-%j.out
#SBATCH --mail-user=g39103001@gm.ym.edu.tw
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
module load gcc/5.4.0
module load viennarna/2.4.9
module load bowtie/1.1.2
module load blast+/2.6.0
#
#= JAVA memory allocation (space for else than RDD operations)
export _JAVA_OPTIONS="-Xms3g -Xmx6g"


## --------------------------------------
## footnotes
## --------------------------------------
#= sbatch pyspark_submit_jv2_02.sh
#= squeue -u cjwu
#= scancel <jobid>
#
#= https://www.sherlock.stanford.edu/docs/software/using/spark/

