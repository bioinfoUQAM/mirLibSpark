'''
program: mirLibPipeline_dummy.py
date: 2018-09-21
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

  if not len(sys.argv) == 2:
    sys.stderr.write('Three arguments required\nUsage: spark-submit mirLibPipeline.py <path to paramfile> 2>/dev/null\n')
    sys.exit()

  paramfile = sys.argv[1]
  paramDict = ut.readParam (paramfile)

  #= Parameters and cutoffs
  input_type = paramDict['input_type']
  adapter = ut.tr_U_T (paramDict['adapter'])
  project_path = paramDict['project_path'][:-1]
  rep_input = paramDict['input_path']
  rep_output = paramDict['output_path']
  rep_msub_jobsOut = project_path + '/workdir/jobsOut'
  rep_tmp = project_path + '/tmp/'                   

  #= spark configuration
  appMaster = paramDict['sc_master']                #"local[*]" 
  appName = paramDict['sc_appname']                 #"mirLibSpark"
  mstrMemory = paramDict['sc_mstrmemory']           #"4g"
  execMemory = paramDict['sc_execmemory']           #"4g"
  execCores = paramDict['sc_execcores']             #2
  partition = int(paramDict['sc_partition'])

  #= genome
  genome_path = paramDict['genome_path'] 

  #= cutoffs
  limit_srna_freq = int(paramDict['limit_s_freq'])  #10       # exclude sRNA freq < limit_srna_freq
  limit_mrna_freq = int(paramDict['limit_m_freq'])  #200      # exclude miRNA freq < limit_mrna_freq
  limit_len = int(paramDict['limit_len'])           #18       # exclude RNA length < limit_len
  limit_nbLoc = int(paramDict['limit_nbLoc'])       #3        # exculde nbLoc mapped with bowtie  > limit_nbLoc

  #= bowtie
  b_index_path = paramDict['b_index_path']
  chromosomes = paramDict['chromosomes'].split(',')
  bowtie_index_suffix = paramDict['bowtie_index_suffix']

  #= file and list of known non miRNA
  known_non = paramDict['known_non_file'] 
  d_ncRNA_CDS = ut.get_nonMirna_coors (known_non) #= nb = 198736

  #= RNAfold
  path_RNAfold = project_path + '/lib/'
  path_RNAfold = ut.find_RNAfold_path () #mirdup dependency
  temperature = int(paramDict['temperature']) 

  #= pri-mirna
  pri_l_flank = int(paramDict['pri_l_flank'])       #120
  pri_r_flank = int(paramDict['pri_r_flank'])       #60
  pre_flank = int(paramDict['pre_flank'])           #30

  #= mircheck parameter
  mcheck_param = paramDict['mcheck_param']          #'def'    # def : default parameters / mey : meyers parameters

  #= miRdup parameter
  mirdup_model = project_path + '/lib/miRdup_1.4/model/' + paramDict['mirdup_model']
  mirdup_jar = project_path + '/lib/miRdup_1.4/miRdup.jar'
  mirdup_limit =  float(paramDict['mirdup_limit'])

  #= miRanda parameter
  target_file = paramDict['target_file']
  miranda_binary = project_path + '/lib/miranda'
  Max_Score_cutoff = paramDict['Max_Score_cutoff'] #= need string or buffer
  Max_Energy_cutoff = paramDict['Max_Energy_cutoff'] #= NOT WORKING YET
  Gap_Penalty = paramDict['Gap_Penalty']

  #= EXMAMINE OPTIONS 
  ut.validate_options(paramDict)

  #= make required folders if not exist
  reps = [rep_output, rep_tmp, rep_msub_jobsOut]
  ut.makedirs_reps (reps)

  #= Objects for rule functions
  dmask_obj = mru.prog_dustmasker()
  dmask_cmd, dmask_env = dmask_obj.dmask_pipe_cmd()
  kn_obj = mru.prog_knownNonMiRNA(d_ncRNA_CDS)
  rnafold_obj = mru.prog_RNAfold(temperature)
  mircheck_obj = mru.prog_mirCheck(mcheck_param)
  mirdup_obj = mru.prog_miRdup (rep_tmp, mirdup_model, mirdup_jar, path_RNAfold)
  profile_obj = mru.prog_dominant_profile()
  miranda_obj = mru.prog_miRanda(Max_Score_cutoff, Max_Energy_cutoff, target_file, rep_tmp, miranda_binary, Gap_Penalty)

  
  #= Spark context
  sc = ut.pyspark_configuration(appMaster, appName, mstrMemory, execMemory, execCores)
  #
  sc.addFile(known_non)
  #
  sc.addPyFile(project_path + '/src/utils.py')
  sc.addPyFile(project_path + '/src/mirLibRules.py')
  sc.addFile(project_path + '/src/eval_mircheck.pl')
  sc.addFile(project_path + '/lib/miRcheck.pm')
  sc.addFile(project_path + '/lib/miRdup_1.4/lib/weka.jar')
  sc.addFile(project_path + '/lib/dustmasker')
  sc.addFile(project_path + '/lib/RNAfold')
  sc.addFile(project_path + '/lib/bowtie')
  sc.addFile(project_path + '/lib/bowtie-align-l')
  sc.addFile(project_path + '/lib/bowtie-align-s')
  sc.addFile(project_path + '/lib/VARNAv3-93.jar')
  sc.addFile(mirdup_jar)
  sc.addFile(mirdup_model)
  sc.addFile(target_file)
  sc.addFile(miranda_binary)

  #= Spark application ID
  appId = str(sc.applicationId)
  
  #= Fetch library files in rep_input
  infiles = [f for f in listdir(rep_input) if os.path.isfile(os.path.join(rep_input, f))]
  
  #= Time processing of libraries
  timeDict = {}
  
  
  for infile in infiles :
    if infile[-1:] == '~': continue
    print ("--Processing of the library: ", infile)

    inBasename = os.path.splitext(infile)[0]
    infile = rep_input+infile

    print ("  Start of the processing...", end="\n")
    startLib = time.time()
    
    distFile_rdd = sc.textFile("file:///" + infile, partition)

    if input_type == 'a': #= raw
      collapse_rdd = distFile_rdd.map(lambda line: mru.rearrange_rule(line, '\t'))

  predictionResults = collapse_rdd.collect() 

  print('predictionResults:')
  print(predictionResults)
  
  sc.stop() #= allow to run multiple SparkContexts



  ### test to initiate a new sc context ########
  sc = ut.pyspark_configuration(appMaster, appName, mstrMemory, execMemory, execCores) 
  distFile_rdd = sc.textFile("file:///" + infile, partition) 
  second_sccontext_results = distFile_rdd.collect() 
  
  print('second_sccontext_results:') 
  print(second_sccontext_results) 
  sc.stop() 

