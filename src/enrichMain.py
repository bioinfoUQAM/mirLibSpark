'''
Chao-Jung Wu
2018-09-25
'''
import os


#cmd = 'git pull origin tgkegg';os.system (cmd)




#========================================
# MAIN
#========================================
#= perl compute_enrichment.pl [namecode.txt] [outputfolder] [option 1 or 2]
#= 1: upper cumulative
#= 2: lower comulative
#= ex: cmd = 'perl compute_enrichment.pl ../output/preXLOC_GO/namecode.txt ../output/preXLOC_GO/output_comput_enrich 1';os.system(cmd)
keyword = 'top5scoredTG_pathways'
#cmd = 'mkdir ' + keyword ;os.system(cmd)
cmd = 'perl compute_enrichment.pl ../output/'+ keyword +'/namecode.txt ../output/'+ keyword +'/output_comput_enrich 1';os.system(cmd)
#========================================





#=========================================
#= cloudera install R and packages
#cmd = 'su';os.system(cmd)
#cmd = 'yum install R';os.system(cmd)
#cmd = 'R';os.system(cmd)
#=within R, at the prompt
#install.packages("RColorBrewer")
#install.packages("gplots")
#q() #= exit R shell
#exit #= exit su
#=========================================




#=========================================
#=https://docs.computecanada.ca/wiki/R
#=Cedra looking for RColorBrewer and gplots
#cmd = 'module spider r/3.5.0';os.system(cmd)
#cmd = 'module load nixpkgs/16.09';os.system(cmd)
#cmd = 'module load gcc/5.4.0';os.system(cmd)
#cmd = 'module load r/3.5.0';os.system(cmd)
#cmd = 'module load circos/0.69-6';os.system(cmd) #== testing if this has RColorBrewer and gplots
#=========================================
'''
#########################################################
# 2018-09-25
# Golrokh
# to convert bam to bed
#########################################################
#!/bin/bash
#SBATCH --job-name=180911_Lim_bam2bed
#SBATCH --account=def-banire
#SBATCH --time=24:30:00
#SBATCH --mem=515000M
#SBATCH --output=/project/def-banire/Labobioinfo/Jobs/180911_Lim/results/bam2bed.out
#SBATCH --error=/project/def-banire/Labobioinfo/Jobs/180911_Lim/results/bam2bed.err
#SBATCH --mail-user=mamoolack@gmail.com
#SBATCH --mail-type=ALL
module load r/3.5.0
#add other modules
#add running the Rscript
'''


