#!/bin/bash
#SBATCH --account=def-banire
#SBATCH --output=jobout/%x-%j.out
#SBATCH --mail-user=wu.chaojung@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --job-name=simple-181012
#SBATCH --time=00:01:00
echo 'Hello, world!'
sleep 30
