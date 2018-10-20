#!/bin/bash
#SBATCH --account=def-banire
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --mem=5000M
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=build-dbs_ensembl40
#SBATCH --mail-user=wu.chaojung@gmail.com
#SBATCH --mail-type=ALL


module load nixpkgs/16.09
module load gcc/5.4.0
module load bowtie/1.1.2

python init_dbs_ensembl40_v2.py wheat 2 build	
python init_dbs_ensembl40_v2.py corn 2 build
python init_dbs_ensembl40_v2.py rice 1 build
python init_dbs_ensembl40_v2.py potato 1 build
python init_dbs_ensembl40_v2.py brome 1 build
python init_dbs_ensembl40_v2.py wheatD 1 build




#= sbatch submit.sh
#= squeue -u cjwu
#= scancel <jobid>


