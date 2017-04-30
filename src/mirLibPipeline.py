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
'''

from __future__ import print_function
import sys
import os.path
import time
from os import listdir
#
import utils as ut
import mirLibRules as mru


if __name__ == '__main__' :

  if not len(sys.argv) == 4:
    sys.stderr.write('Three arguments required\nUsage: spark-submit mirLibPipeline.py <path to paramfile> <path to the input directory> <path to the output directory>2>/dev/null\n')
    sys.exit()

  paramfile = sys.argv[1]
  mypath = sys.argv[2]
  rep_output = sys.argv[3]
  
  paramDict = ut.readParam (paramfile)

  # Parameters and cutoffs
  my_sep = paramDict['my_sep']                      # Separator
  rep_tmp = paramDict['rep_tmp']                    # tmp file folder
  # spark configuration
  appMaster = paramDict['sc_master']                #"local" 
  appName = paramDict['sc_appname']                 #"mirLibHadoop"
  mstrMemory = paramDict['sc_mstrmemory']           #"4g"
  execMemory = paramDict['sc_execmemory']           #"4g"
  execNb = paramDict['sc_execnb']                   #4
  execCores = paramDict['sc_execcores']             #2
  # genome
  genome_path = paramDict['genome_path']            #"../input/ATH/TAIR/Genome/"
  # cutoffs
  limit_srna_freq = int(paramDict['limit_s_freq'])  #10       # exclude sRNA freq < limit_srna_freq
  limit_mrna_freq = int(paramDict['limit_m_freq'])  #200      # exclude miRNA freq < limit_mrna_freq
  limit_len = int(paramDict['limit_len'])           #18       # exclude RNA length < limit_len
  limit_nbLoc = int(paramDict['limit_nbLoc'])       #3        # exculde nbLoc mapped with bowtie  > limit_nbLoc
  # bowtie
  b_index = paramDict['b_index']
  # pri-mirna
  pri_l_flank = int(paramDict['pri_l_flank'])       #120
  pri_r_flank = int(paramDict['pri_r_flank'])       #60
  pre_flank = int(paramDict['pre_flank'])           #30
  # mircheck parameter
  mcheck_param = paramDict['mcheck_param']          #'def'    # def : default parameters / mey : meyers parameters

  # Spark context
  sc = ut.pyspark_configuration(appMaster, appName, mstrMemory, execMemory, execNb, execCores)
  sc.addPyFile('utils.py')
  sc.addPyFile('mirLibRules.py')
  # Spark application ID
  appId = str(sc.applicationId)
  
  # Objects for rule functions
  dmask_obj = mru.prog_dustmasker()
  dmask_cmd, dmask_env = dmask_obj.dmask_pipe_cmd()
  bowtie_obj = mru.prog_bowtie(b_index)
  bowtie_cmd, bowtie_env = bowtie_obj.Bowtie_pipe_cmd()
  prec_obj = mru.extract_precurosrs(genome_path, pri_l_flank, pri_r_flank, pre_flank)
  rnafold_obj = mru.prog_RNAfold()
  mircheck_obj = mru.prog_mirCheck(mcheck_param)
  profile_obj = mru.prog_dominant_profile()

  # Fetch library files in mypath
  infiles = [f for f in listdir(mypath) if os.path.isfile(os.path.join(mypath, f))]
  
  # Time processing of libraries
  timeDict = {}
  
  for infile in infiles :
    print ("--Processing of the library: ", infile)
    
    infile = mypath+infile
    inBasename = os.path.splitext(os.path.basename(infile))[0]
    inKvfile = rep_tmp + inBasename + '.kv.txt'
    hdfsFile = inBasename + '.hkv.txt'

    # Convert the input file to a Key value file
    ut.convert_seq_freq_file_to_KeyValue(infile, inKvfile, my_sep)
  
    # Save a local file to HDFS system
    ut.convertTOhadoop(inKvfile, hdfsFile)
    
    #
    print ("  Start of the processing...", end="\n")
    startLib = time.time()
    # Convert the text file to RDD object
    distFile = sc.textFile(hdfsFile)
    input_rdd = distFile.map(lambda line: mru.rearrange_rule(line, my_sep)) # (seq, freq)
  
    # Filtering sRNA low frequency
    sr_low_rdd = input_rdd.filter(lambda e: int(e[1]) > limit_srna_freq)

    # Filtering short length
    sr_short_rdd = sr_low_rdd.filter(lambda e: len(e[0]) > limit_len).persist()

    # Filtering with DustMasker
    dmask_rdd = sr_short_rdd.map(lambda e: '>s\n'+e[0])\
                            .pipe(dmask_cmd, dmask_env)\
                            .filter(lambda e: e.isupper() and not e.startswith('>'))\
                            .map(lambda e: str(e.rstrip()))\
                            .persist()
    
    # Mapping with Bowtie
    bowtie_rdd = dmask_rdd.pipe(bowtie_cmd, bowtie_env)\
                          .map(bowtie_obj.bowtie_rearrange_map)\
                          .groupByKey()\
                          .map(lambda e: (e[0], [len(list(e[1])), list(e[1])]))\
                          .persist()
    
    # Get the expression value for each reads
    bowFrq_rdd = bowtie_rdd.join(sr_short_rdd)\
                           .map(bowtie_obj.bowtie_freq_rearrange_rule)\
                           .persist()
    
    # Filtering miRNA low frequency
    mr_low_rdd = bowFrq_rdd.filter(lambda e: e[1][0] > limit_mrna_freq)
    
    # Filtering high nbLocations and zero location
    nbLoc_rdd = mr_low_rdd.filter(lambda e: e[1][1] > 0 and e[1][1] < limit_nbLoc)
    
    # Extraction of the pri-miRNA
    primir_rdd = nbLoc_rdd.flatMap(prec_obj.extract_prim_rule)
    
    # pri-miRNA folding
    pri_fold_rdd = primir_rdd.map(lambda e: rnafold_obj.RNAfold_map_rule(e, 3))
    
    # Validating pri-mirna with mircheck
    pri_vld_rdd = pri_fold_rdd.map(lambda e: mircheck_obj.mirCheck_map_rule(e, 3))\
                              .filter(lambda e: any(e[1][3]))
    
    # Filtering structure with branched loop
    one_loop_rdd = pri_vld_rdd.filter(lambda e: ut.containsOnlyOneLoop (e[1][3][2][int(e[1][3][4]) : int(e[1][3][5])+1]))
    
    # Extraction of the pre-miRNA
    premir_rdd = one_loop_rdd.map(lambda e: prec_obj.extract_prem_rule(e, 3))

    # pre-miRNA folding
    pre_fold_rdd = premir_rdd.map(lambda e: rnafold_obj.RNAfold_map_rule(e, 4))
    
    # Validating pre-mirna with mircheck
    pre_vld_rdd = pre_fold_rdd.map(lambda e: mircheck_obj.mirCheck_map_rule(e, 4))\
                              .filter(lambda e: any(e[1][4]))
    
    # you can use chromo_strand as key to search bowtie blocs in the following dict
    dict_bowtie_chromo_strand = profile_obj.get_bowtie_strandchromo_dict(bowFrq_rdd.collect())
    
    # Results of miRNA prediction
    #miRNA_rdd = pre_vld_rdd.filter(lambda e: profile_obj.exp_profile_filter(e, dict_bowtie_chromo_strand))

    miRNA_rdd = pre_vld_rdd.map(lambda e: profile_obj.sudo(e, dict_bowtie_chromo_strand))\
                      .filter(lambda e: e[1][0] / float(e[1][5]) > 0.2)

    #'''
    results = miRNA_rdd.collect()
    
    #
    endLib = time.time()
    print ("  End of the processing     ", end="\n")
    
    # write results to a file
    outFile = rep_output + inBasename + '_miRNAprediction.txt'
    ut.writeToFile (results, outFile)
    #'''
    
    timeDict[inBasename] = endLib - startLib
    
  sc.stop() #allow to run multiple SparkContexts
  
  # print executions time  to a file
  outTime = rep_output + appId + '_time.txt'
  ut.writeTimeLibToFile (timeDict, outTime, appId, paramDict)
