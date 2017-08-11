'''
program: mirLibPipeline.py
author: M.A.Remita
author: Chao-Jung Wu
date: 2017-03-28
version: 1.00.02

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

  if not len(sys.argv) == 2:
    sys.stderr.write('Three arguments required\nUsage: spark-submit mirLibPipeline.py <path to paramfile> 2>/dev/null\n')
    sys.exit()

  paramfile = sys.argv[1]
  paramDict = ut.readParam (paramfile)

  #= Parameters and cutoffs
  #platform = paramDict['platform'] 
  #project_path = paramDict['project_path_' + platform][:-1]
  project_path = paramDict['project_path'][:-1]
  rep_input = paramDict['input_path']
  rep_output = paramDict['output_path']
  rep_msub_jobsOut = project_path + '/workdir/jobsOut'
  my_sep = paramDict['my_sep']                      # Separator
  rep_tmp = project_path + '/tmp/'                   # tmp file folder
  #= spark configuration

  appMaster = paramDict['sc_master']                #"local" 
  appName = paramDict['sc_appname']                 #"mirLibHadoop"
  mstrMemory = paramDict['sc_mstrmemory']           #"4g"
  execMemory = paramDict['sc_execmemory']           #"4g"
  execNb = paramDict['sc_execnb']                   #4
  execCores = paramDict['sc_execcores']             #2

  #= genome
  #genome_path = paramDict['genome_path_' + platform] 
  genome_path = paramDict['genome_path'] 
  #= cutoffs

  limit_srna_freq = int(paramDict['limit_s_freq'])  #10       # exclude sRNA freq < limit_srna_freq
  limit_mrna_freq = int(paramDict['limit_m_freq'])  #200      # exclude miRNA freq < limit_mrna_freq
  limit_len = int(paramDict['limit_len'])           #18       # exclude RNA length < limit_len
  limit_nbLoc = int(paramDict['limit_nbLoc'])       #3        # exculde nbLoc mapped with bowtie  > limit_nbLoc

  # bowtie
  b_index_path = paramDict['b_index_path']

  #= file and list of known non miRNA
  known_non = '../dbs/TAIR10_ncRNA_CDS.gff' ###############################################################################
  #l_non_miRNA = ut.get_nonMirna_list (known_non, genome_path)
  #l_non_miRNA = ['TGGATTTATGAAAGACGAACAACTGCGAAA']
  d_ncRNA_CDS = ut.get_nonMirna_coors (known_non) #= nb = 198736
  #d_ncRNA_CDS = {1: ['+', 'Chr5', '26939753', '26939884'], 2: ['+', 'Chr5', '26939972', '26940240'], 3: ['+', 'Chr5', '26940312', '26940396'], 4: ['+', 'Chr5', '26940532', '26940578']}

  # pri-mirna
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
  #= miRanda parameter
  #target_file = project_path + '/lib/' + paramDict['target_file']
  #miranda_exe = project_path + '/lib/miranda'
  #Max_Score_cutoff = float(paramDict['Max_Score_cutoff'])
  #query_motif_match_cutoff = float(paramDict['query_motif_match_cutoff'])
  #gene_motif_match_cutoff = float(paramDict['gene_motif_match_cutoff'])
  #Max_Energy_cutoff = float(paramDict['Max_Energy_cutoff'])
  #= make required folders if not exist
  reps = [rep_output, rep_tmp, rep_msub_jobsOut]
  ut.makedirs_reps (reps)
  
  #= Spark context

  sc = ut.pyspark_configuration(appMaster, appName, mstrMemory, execMemory, execNb, execCores)
  #
  sc.addPyFile(project_path + '/src/utils.py')
  sc.addPyFile(project_path + '/src/mirLibRules.py')
  sc.addFile(project_path + '/src/eval_mircheck.pl')
  sc.addFile(project_path + '/lib/miRcheck.pm')

  sc.addFile(project_path + '/lib/miRdup_1.4/lib/weka.jar')
  sc.addFile(mirdup_jar)
  sc.addFile(mirdup_model)

  #sc.addFile(target_file)
  #sc.addFile(miranda_exe)

  #= Spark application ID

  appId = str(sc.applicationId)
  
  #= Objects for rule functions
  dmask_obj = mru.prog_dustmasker()
  dmask_cmd, dmask_env = dmask_obj.dmask_pipe_cmd()
  bowtie_obj = mru.prog_bowtie(b_index_path)
  kn_obj = mru.prog_knownNonMiRNA(d_ncRNA_CDS) ##########################################################
  bowtie_cmd, bowtie_env = bowtie_obj.Bowtie_pipe_cmd()
  prec_obj = mru.extract_precurosrs(genome_path, pri_l_flank, pri_r_flank, pre_flank)
  rnafold_obj = mru.prog_RNAfold()
  mircheck_obj = mru.prog_mirCheck(mcheck_param)
  profile_obj = mru.prog_dominant_profile()

  mirdup_obj = mru.prog_miRdup (rep_tmp, mirdup_model, mirdup_jar, path_RNAfold)
  #miranda_obj = mru.prog_miRanda(Max_Score_cutoff, query_motif_match_cutoff, gene_motif_match_cutoff, Max_Energy_cutoff, target_file, rep_tmp, miranda_exe)

  #= Fetch library files in rep_input

  infiles = [f for f in listdir(rep_input) if os.path.isfile(os.path.join(rep_input, f))]
  
  #= Time processing of libraries
  timeDict = {}
  
  for infile in infiles :
    if infile[:-1] == '~': continue
    print ("--Processing of the library: ", infile)
    
    infile = rep_input+infile
    inBasename = os.path.splitext(os.path.basename(infile))[0]
    inKvfile = rep_tmp + inBasename + '.kv.txt'

    # hdfsFile = inBasename + '.hkv.txt'

    #= Convert the input file to a Key value file
    ut.convert_seq_freq_file_to_KeyValue(infile, inKvfile, my_sep)
    
    #
    print ("  Start of the processing...", end="\n")
    startLib = time.time()
    
    #= Convert the text file to RDD object
    ## in : file
    ## out: u'seq,freq'
    distFile = sc.textFile("file:///" + inKvfile)
   
    #= Convert element from string to list
    ## in : u'seq,freq'
    ## out: ('seq', freq)
    input_rdd = distFile.map(lambda line: mru.rearrange_rule(line, my_sep)).persist()
    print('NB input_rdd: ', len(input_rdd.collect()))####################################################
   
    #= Filtering sRNA low frequency
    ## in : ('seq', freq)
    ## out: ('seq', freq)
    sr_low_rdd = input_rdd.filter(lambda e: int(e[1]) > limit_srna_freq).persist()#########################
    print('NB sr_low_rdd: ', len(sr_low_rdd.collect()))####################################################
    
    #= Filtering short length
    ## in : ('seq', freq)
    ## out: ('seq', freq)
    sr_short_rdd = sr_low_rdd.filter(lambda e: len(e[0]) > limit_len).persist()  # TO KEEP IT
    print('NB sr_short_rdd: ', len(sr_short_rdd.collect()))###################################################
    
    #= Filtering with DustMasker
    ## in : ('seq', freq)
    ## out: 'seq'
    dmask_rdd = sr_short_rdd.map(lambda e: '>s\n'+e[0])\
                            .pipe(dmask_cmd, dmask_env)\
                            .filter(lambda e: e.isupper() and not e.startswith('>'))\
                            .map(lambda e: str(e.rstrip()))\
                            .persist()
    print('NB dmask_rdd: ', len(dmask_rdd.collect()))############################################
    
    #= Mapping with Bowtie
    ## in : 'seq'
    ## out: ('seq', [nbLoc, [['strd','chr',posChr],..]])
    bowtie_rdd = dmask_rdd.map(lambda e: " "+e)\
                          .pipe(bowtie_cmd, bowtie_env)\
                          .map(bowtie_obj.bowtie_rearrange_map)\
                          .groupByKey()\
                          .map(lambda e: (e[0], [len(list(e[1])), list(e[1])]))\
                          .persist()
    print('NB bowtie_rdd: ', len(bowtie_rdd.collect()))##################################################
   
    #= Getting the expression value for each reads
    ## in : ('seq', [nbLoc, [['strd','chr',posChr],..]])
    ## out: ('seq', [freq, nbLoc, [['strd','chr',posChr],..]])
    bowFrq_rdd = bowtie_rdd.join(sr_short_rdd)\
                           .map(bowtie_obj.bowtie_freq_rearrange_rule)\
                           .persist()
    print('NB bowFrq_rdd: ', len(bowFrq_rdd.collect()))##################################################
    
    #= Filtering miRNA low frequency
    ## in : ('seq', [freq, nbLoc, [['strd','chr',posChr],..]])
    ## out: ('seq', [freq, nbLoc, [['strd','chr',posChr],..]])
    mr_low_rdd = bowFrq_rdd.filter(lambda e: e[1][0] > limit_mrna_freq).persist()########################
    print('NB mr_low_rdd: ', len(mr_low_rdd.collect()))##################################################
    
    #= Filtering high nbLocations and zero location
    ## in : ('seq', [freq, nbLoc, [['strd','chr',posChr],..]])
    ## out: ('seq', [freq, nbLoc, [['strd','chr',posChr],..]])
    nbLoc_rdd = mr_low_rdd.filter(lambda e: e[1][1] > 0 and e[1][1] < limit_nbLoc).persist()#############
    print('NB nbLoc_rdd: ', len(nbLoc_rdd.collect()))###################################################
    
    #= Flatmap the RDD
    ## in : ('seq', [freq, nbLoc, [['strd','chr',posChr],..]])
    ## out: ('seq', [freq, nbLoc, ['strd','chr',posChr])
    flat_rdd = nbLoc_rdd.flatMap(mru.flatmap_mappings).persist() ###
    print('NB flat_rdd (this step flats elements): ', len(flat_rdd.groupByKey().collect()))##################

    ###############################
    ## Filtering known non-miRNA ##
    ###############################
    ## in : ('seq', [freq, nbLoc, ['strd','chr',posChr])
    ## out: ('seq', [freq, nbLoc, ['strd','chr',posChr])
    ##excluKnownNon_rdd = flat_rdd.filter(kn_obj.knFilterBySeq) #= defunct
    excluKnownNon_rdd = flat_rdd.filter(kn_obj.knFilterByCoor).persist()#######
    print('excluKnownNon_rdd: ', len(excluKnownNon_rdd.groupByKey().collect()))########
    
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
    print('NB pri_vld_rdd distinct (mircheck): ', len(pri_vld_rdd.groupByKey().collect()))#################################

    #= Filtering structure with branched loop
    ## in : ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold','mkPred','mkStart','mkStop']])
    ## out: ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold','mkPred','mkStart','mkStop']])
    one_loop_rdd = pri_vld_rdd.filter(lambda e: ut.containsOnlyOneLoop(e[1][3][2][int(e[1][3][4]) : int(e[1][3][5])+1])).persist()############
    print('NB one_loop_rdd distinct : ', len(one_loop_rdd.groupByKey().collect()))#########################
    
    #= Extraction of the pre-miRNA
    ## in : ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold','mkPred','mkStart','mkStop']])
    ## out: ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold','mkPred','mkStart','mkStop'], ['preSeq',posMirPre]])
    premir_rdd = one_loop_rdd.map(lambda e: prec_obj.extract_prem_rule(e, 3))
    
    #= pre-miRNA folding
    ## in : ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold', 'mkPred','mkStart','mkStop'], ['preSeq',posMirPre]])
    ## out: ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold', 'mkPred','mkStart','mkStop'], ['preSeq',posMirPre,'preFold']])
    pre_fold_rdd = premir_rdd.map(lambda e: rnafold_obj.RNAfold_map_rule(e, 4))

    ###################################################   
    #= Validating pre-mirna with mircheck
    ## in : ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold', 'mkPred','mkStart','mkStop'], ['preSeq',posMirPre,'preFold']])
    ## out: ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold', 'mkPred','mkStart','mkStop'], ['preSeq',posMirPre,'preFold','mkPred','mkStart','mkStop']])
    pre_vld_rdd0 = pre_fold_rdd.map(lambda e: mircheck_obj.mirCheck_map_rule(e, 4))\
                              .filter(lambda e: any(e[1][4])).persist()
    print('NB pre_vld_rdd0 distinct (mircheck II): (', len(pre_vld_rdd0.groupByKey().collect()), ')')################

    #= Validating pre-mirna with miRdup zipWithUniqueId
    ## in : ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold', 'mkPred','mkStart','mkStop'], ['preSeq',posMirPre,'preFold']])
    ## out: ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold', 'mkPred','mkStart','mkStop'], ['preSeq',posMirPre,'preFold','mpPred','mpScore']])
    pre_vld_rdd = pre_fold_rdd.zipWithIndex()\
                              .map(mirdup_obj.run_miRdup)\
                              .map(lambda e: e[0])\
                              .filter(lambda e: e[1][4][3] == "true")\
                              .persist()##################
    print('NB pre_vld_rdd distinct (mirdup): ', len(pre_vld_rdd.groupByKey().collect()))################
    ###################################################

    #= Create dict, chromo_strand as key to search bowtie blocs in the following dict
    dict_bowtie_chromo_strand = profile_obj.get_bowtie_strandchromo_dict(bowFrq_rdd.collect())
    
    #= Filtering by expression profile (< 20%)
    ## in : ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold', 'mkPred','mkStart','mkStop'], ['preSeq',posMirPre,'preFold','mpPred','mpScore']])
    ## out: ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold', 'mkPred','mkStart','mkStop'], ['preSeq',posMirPre,'preFold','mpPred','mpScore'], totalfrq])
    profile_rdd = pre_vld_rdd.map(lambda e: profile_obj.computeProfileFrq(e, dict_bowtie_chromo_strand))\
                      .filter(lambda e: e[1][0] / float(e[1][5]) > 0.2)\
                      .persist()####################
    print('NB profile_rdd distinct (final prediction): ', len(profile_rdd.groupByKey().collect()))#####################

    #= target prediction
    #miranda_rdd = miRNA_rdd.map(miranda_obj.computeTargetbyMiranda).persist()####
    #print('NB miranda_rdd distinct : ', len(miranda_rdd.groupByKey().collect()))####
    ##results = miranda_rdd.collect()

    #results = miRNA_rdd.collect()

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
