#!/bin/bash
#SBATCH --job-name=mirdeepp
#SBATCH --account=youraccount
#SBATCH --mail-user=yourname@email.com
#SBATCH --time=60:10:00
#SBATCH --nodes=1
#SBATCH --mem=5000M
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --error=%x-%j.err
#SBATCH --output=%x-%j.out
#SBATCH --mail-type=ALL

module load nixpkgs/16.09
module load gcc/5.4.0
module load viennarna/2.4.9
module load bowtie/1.1.2
module load blast+/2.6.0
module load python/2.7.14

time python run_mirdeep_p.py 2>/dev/null
