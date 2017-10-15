'''
program: mirLibPipeline_wheat.py
author: M.A.Remita
author: Chao-Jung Wu
date: 2017-10-14
version: 1.00.00

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
    sys.stderr.write('Three arguments required\nUsage: spark-submit mirLibPipeline_wheat.py <path to paramfile> 2>/dev/null\n')
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
  #my_sep = paramDict['my_sep']                      # Separator
  rep_tmp = project_path + '/tmp/'                   # tmp file folder

  #= spark configuration
  appMaster = paramDict['sc_master']                #"local" 
  appName = paramDict['sc_appname']                 #"mirLibHadoop"
  mstrMemory = paramDict['sc_mstrmemory']           #"4g"
  execMemory = paramDict['sc_execmemory']           #"4g"
  execNb = paramDict['sc_execnb']                   #4
  execCores = paramDict['sc_execcores']             #2

  #= genome
  genome_path = paramDict['genome_path'] 

  #= cutoffs
  limit_srna_freq = int(paramDict['limit_s_freq'])  #10       # exclude sRNA freq < limit_srna_freq
  limit_mrna_freq = int(paramDict['limit_m_freq'])  #200      # exclude miRNA freq < limit_mrna_freq
  limit_len = int(paramDict['limit_len'])           #18       # exclude RNA length < limit_len
  limit_nbLoc = int(paramDict['limit_nbLoc'])       #3        # exculde nbLoc mapped with bowtie  > limit_nbLoc

  #= bowtie
  b_index_path = paramDict['b_index_path']

  #= file and list of known non miRNA
  #known_non = '../dbs/TAIR10_ncRNA_CDS.gff'
  known_non = project_path + '/dbs/TAIR10_ncRNA_CDS.gff'###########################
  d_ncRNA_CDS = ut.get_nonMirna_coors (known_non) #= nb = 198736

  #= pri-mirna
  pri_l_flank = int(paramDict['pri_l_flank'])       #120
  pri_r_flank = int(paramDict['pri_r_flank'])       #60
  pre_flank = int(paramDict['pre_flank'])           #30

  #= mircheck parameter
  mcheck_param = paramDict['mcheck_param']          #'def'    # def : default parameters / mey : meyers parameters

  #= miRdup parameter
  path_RNAfold = ut.find_RNAfold_path ()
  mirdup_model = project_path + '/lib/miRdup_1.4/model/' + paramDict['mirdup_model']
  mirdup_jar = project_path + '/lib/miRdup_1.4/miRdup.jar'
  mirdup_limit =  float(paramDict['mirdup_limit'])

  ## EXMAMIN OPTIONS ####################################
  ut.validate_options(paramDict)
  #######################################################

  #= make required folders if not exist
  reps = [rep_output, rep_tmp, rep_msub_jobsOut]
  ut.makedirs_reps (reps)
  
  #= Spark context
  sc = ut.pyspark_configuration(appMaster, appName, mstrMemory, execMemory, execNb, execCores)
  #
  sc.addFile(known_non)################################
  sc.addPyFile(project_path + '/src/utils.py')
  sc.addPyFile(project_path + '/src/mirLibRules.py')
  sc.addFile(project_path + '/src/eval_mircheck.pl')
  sc.addFile(project_path + '/lib/miRcheck.pm')
  sc.addFile(project_path + '/lib/miRdup_1.4/lib/weka.jar')
  sc.addFile(mirdup_jar)
  sc.addFile(mirdup_model)

  #= Spark application ID
  appId = str(sc.applicationId)
  
  #= Objects for rule functions
  dmask_obj = mru.prog_dustmasker()
  dmask_cmd, dmask_env = dmask_obj.dmask_pipe_cmd()
  bowtie_obj = mru.prog_bowtie(b_index_path)
  kn_obj = mru.prog_knownNonMiRNA(d_ncRNA_CDS)
  bowtie_cmd, bowtie_env = bowtie_obj.Bowtie_pipe_cmd()
  prec_obj = mru.extract_precurosrs(genome_path, pri_l_flank, pri_r_flank, pre_flank)
  rnafold_obj = mru.prog_RNAfold()
  mircheck_obj = mru.prog_mirCheck(mcheck_param)
  profile_obj = mru.prog_dominant_profile()

  mirdup_obj = mru.prog_miRdup (rep_tmp, mirdup_model, mirdup_jar, path_RNAfold)

  #= Fetch library files in rep_input
  infiles = [f for f in listdir(rep_input) if os.path.isfile(os.path.join(rep_input, f))]
  
  #= Time processing of libraries
  timeDict = {}
  
  for infile in infiles :
    if infile[-1:] == '~': continue
    print ("--Processing of the library: ", infile)
    
    inBasename = os.path.splitext(infile)[0]
    infile = rep_input+infile
    inKvfile = rep_tmp + inBasename + '.kv.txt'

    if input_type == 'd': #= fastq
      ut.convert_fastq_file_to_KeyValue(infile, inKvfile)
      infile = inKvfile
      
    #
    print ("  Start of the processing...", end="\n")
    startLib = time.time()
    
    #= Convert the text file to RDD object
    ## in : file
    ## out: (a) u'seq\tfreq', 
    ##      (b) u'seq1', u'seq2', u'seq1', 
    ##      (c) u'>name1\nseq1', u'>name2\nseq2', u'>name3\nseq1',
    ##      (d) u'seq\tquality'
    ##      (e) SAM file
    distFile_rdd = sc.textFile("file:///" + infile, 2) #= NumPartitions = 4 (default for 100.txt was 2)
    #print('distFile_rdd nbPartition', distFile_rdd.getNumPartitions())
    #print('NB distFile_rdd: ', len(distFile_rdd.collect()))#

    if input_type == 'e': #= sam
      sam_rdd = distFile_rdd.filter(lambda line: not line[0] == '@')\
                            .map(lambda line: mru.rearrange_sam_rule(line))\
                            .filter(lambda e: not e[1][2][1][0] == '*')


      excluKnownNon_rdd = sam_rdd
                                
    #= Extraction of the pri-miRNA
    ## in : ('seq', [freq, nbLoc, ['strd','chr',posChr])
    ## out: ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri]])
    primir_rdd = excluKnownNon_rdd.flatMap(prec_obj.extract_prim_rule)

    #= pri-miRNA folding
    ## in : ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri]])
    ## out: ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold']])
    pri_fold_rdd = primir_rdd.map(lambda e: rnafold_obj.RNAfold_map_rule(e, 3))
    
    #= Validating pri-mirna with mircheck
    ## in : ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold']])
    ## out: ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold','mkPred','mkStart','mkStop']])
    pri_vld_rdd = pri_fold_rdd.map(lambda e: mircheck_obj.mirCheck_map_rule(e, 3))\
                              .filter(lambda e: any(e[1][3])).persist()###################
    #print('NB pri_vld_rdd distinct (mircheck): ', len(pri_vld_rdd.groupByKey().collect()))#################################

    #= Filtering structure with branched loop
    ## in : ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold','mkPred','mkStart','mkStop']])
    ## out: ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold','mkPred','mkStart','mkStop']])
    one_loop_rdd = pri_vld_rdd.filter(lambda e: ut.containsOnlyOneLoop(e[1][3][2][int(e[1][3][4]) : int(e[1][3][5])+1]))#.persist()############
    #print('NB one_loop_rdd distinct : ', len(one_loop_rdd.groupByKey().collect()))#########################
    
    #= Extraction of the pre-miRNA
    ## in : ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold','mkPred','mkStart','mkStop']])
    ## out: ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold','mkPred','mkStart','mkStop'], ['preSeq',posMirPre]])
    premir_rdd = one_loop_rdd.map(lambda e: prec_obj.extract_prem_rule(e, 3))
    
    #= pre-miRNA folding
    ## in : ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold', 'mkPred','mkStart','mkStop'], ['preSeq',posMirPre]])
    ## out: ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold', 'mkPred','mkStart','mkStop'], ['preSeq',posMirPre,'preFold']])
    pre_fold_rdd = premir_rdd.map(lambda e: rnafold_obj.RNAfold_map_rule(e, 4))

    #= Validating pre-mirna with miRdup zipWithUniqueId
    ## in : ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold', 'mkPred','mkStart','mkStop'], ['preSeq',posMirPre,'preFold']])
    ## out: ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold', 'mkPred','mkStart','mkStop'], ['preSeq',posMirPre,'preFold','mpPred','mpScore']])
    pre_vld_rdd = pre_fold_rdd.zipWithIndex()\
                              .map(mirdup_obj.run_miRdup)\
                              .filter(lambda e: e[1][4][3] == "true")\
                              .persist()##################
    #print('pre_vld_rdd nbPartition', pre_vld_rdd.getNumPartitions())
    print(pre_vld_rdd.collect())

    
    #'''
    endLib = time.time()
    print ("  End of the processing     ", end="\n")
    
    #= write results to a file
    outFile = rep_output + inBasename + '_miRNAprediction.txt'
    #ut.writeToFile (results, outFile)
    
    
    timeDict[inBasename] = endLib - startLib
    
  sc.stop() #= allow to run multiple SparkContexts
  
  #= print executions time  to a file
  outTime = rep_output + appId + '_time.txt'
  ut.writeTimeLibToFile (timeDict, outTime, appId, paramDict)
