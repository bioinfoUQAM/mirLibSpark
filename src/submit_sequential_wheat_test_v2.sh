#!/bin/bash
#SBATCH --job-name=sequential-181028
#SBATCH --account=def-banire
#SBATCH --time=20:30:00
#SBATCH --nodes=1
#SBATCH --mem=64000M
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --error=jobout_sequential/%x-%j.err
#SBATCH --output=jobout_sequential/%x-%j.out
#SBATCH --mail-user=g39103001@gm.ym.edu.tw
#SBATCH --mail-type=ALL


module load nixpkgs/16.09
module load spark/2.3.0
module load gcc/5.4.0
module load viennarna/2.4.9
module load bowtie/1.1.2
module load blast+/2.6.0


time python sequential_wheat.py

