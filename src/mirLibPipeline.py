'''
program: mirLibPipeline.py
author: M.A.Remita
author: Chao-Jung Wu
date: 2017-03-28
version: 0.00.01

Le programme implemente le pipeline d'analyse des sRAN et prediction des miRNAs. Les principales etapes de predictions sont :
  - 1 Filtrage
  - 2 Alignement 
  - 3 Extraction du precurseur
  - 4 Repliement du precurseur
  - 5 prediction/validation du miRNA
  - 6 validaiton/filtrage avec l'expression
  
La version actuelle accepte un seul argument qui est le fichier contenant les sequences reads a traiter.

personal note - julie's launching command line in her computer:

spark-submit mirLibPipeline.py a_thaliana /home/cloudera/workspace/miRNA_predictor/sudoData/a_th_3.txt 2>/dev/null

spark-submit mirLibPipeline.py a_thaliana /home/cloudera/workspace/miRNA_predictor/sudoData/a_th_10.txt 2>/dev/null

$ time spark-submit mirLibPipeline.py a_thaliana /home/cloudera/workspace/miRNA_predictor/sudoData/110.txt 1>170406_1result_pipeline.txt 2>/dev/null
around 10 mins (6m, 12m, 9.5min)

$ time spark-submit mirLibPipeline.py a_thaliana /home/cloudera/workspace/miRNA_predictor/sudoData/1102.txt 1>/home/cloudera/workspace/miRNA_predictor/logOutput/170407_result_pipeline_1.txt 2>/dev/null

##########################
spark-submit mirLibPipeline.py ../paramfile_julie.txt /home/cloudera/workspace/miRNA_predictor/sudoData/a_th_3.txt 2>/dev/null
'''

import sys

#
import utils as ut
import mirLibRules as mru
import mirTaskMain as task



from os import listdir
from os.path import isfile, join







if __name__ == '__main__' :

  if not len(sys.argv) == 3:
    sys.stderr.write('Two arguments required\nUsage: spark-submit mirLibPipeline.py <path to paramfile> <path to your input> 2>/dev/null\n')
    sys.exit()

  paramfile = sys.argv[1]
  paramDict = ut.readparam (paramfile)

  # batch option
  batch = paramDict['batch']
  if not batch == 'true':
    print 'launching one file'
    infile = sys.argv[2]
    task.SparkTask(paramDict, infile)
  else:
    print 'launching batch files'
    mypath = sys.argv[2]
    infiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    for infile in infiles:
      task.SparkTask(paramDict, mypath+infile)
    



