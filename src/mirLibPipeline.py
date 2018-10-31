'''
program: mirLibPipeline.py
author: M.A.Remita
author: Chao-Jung Wu
date: 2017-03-28
version: 1.01.00

update: 2018-09-20 cjwu, include varna, miranda, html, format output, reduce dependency of third party programs

Le programme implemente le pipeline d'analyse des sRAN et prediction des miRNAs. Les principales etapes de predictions sont :
  - 1 Filtrage
  - 2 Alignement 
  - 3 Extraction du precurseur
  - 4 Repliement du precurseur
  - 5 prediction/validation du miRNA
  - 6 validaiton/filtrage avec l'expression

#= programming note 181028:
#= repartition(num) is better than partitionBy(num) because the latter we might consider segregating chromosomes. But the number of chrs are different in each species. So it is not easyo to say a number that suits everyone.

'''

from __future__ import print_function
import sys
import os.path
import time
import datetime
from os import listdir
#
import utils as ut
import mirLibRules as mru
import arg

#= 1: display intermediate rdd.count(), this makes the time longer
#= 0: not reporting rdd.count() makes the time shorter
reporting = 0 

if __name__ == '__main__' :

  print('\nInitiating and verifying parameters ...')
  paramDict = arg.init_params ()
  print('============================================================\n')
  for k, v in sorted(paramDict.items()): print(k, ': ', v)

  #= spark configuration
  appMaster = paramDict['sc_master']                #"local[*]" 
  appName = paramDict['sc_appname']                 #"mirLibSpark"
  mstrMemory = paramDict['sc_mstrmemory']           #"4g"
  execMemory = paramDict['sc_execmemory']           #"4g"
  #execCores = paramDict['sc_execcores']             #2
  partition = int(paramDict['sc_partition'])

  #= Spark context
  sc = ut.pyspark_configuration(appMaster, appName, mstrMemory, execMemory)

  #= Spark application ID
  appId = str(sc.applicationId)
  #print('spark.executor.memory: ', sc._conf.get('spark.executor.memory'))
  #print('spark.driver.memory: ', sc._conf.get('spark.driver.memory'))
  #print('spark.master: ', sc._conf.get('spark.master'))
  #print('spark.driver.memoryOverhead: ', sc._conf.get('spark.driver.memoryOverhead')) = none
  #print('spark.executor.memoryOverhead: ', sc._conf.get('spark.executor.memoryOverhead')) = none
  #print('spark.cores.max: ', sc._conf.get('spark.cores.max'))

  #= broadcast paramDict
  broadcastVar_paramDict = sc.broadcast(paramDict)
  paramDict = broadcastVar_paramDict.value

  #= Parameters and cutoffs =========================
  #= paths
  input_type = paramDict['input_type']
  adapter = ut.tr_U_T (paramDict['adapter'])
  project_path = paramDict['project_path'][:-1]
  rep_input = paramDict['input_path']
  rep_output = paramDict['output_path'] + '/' + appId + '/'
  #rep_tmp = '../tmp' + appId + '/' #= do not write files in worker node, often mirdup has problem, OUT-OF-MEMORY
  rep_tmp = project_path + '/tmp' + '/' + appId + '/'

  #= print appId to a file
  #outfile = project_path + '/appId.txt'
  #with open (outfile, 'w') as fh: print(appId, file=fh) 

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
  broadcastVar_d_ncRNA_CDS = sc.broadcast(d_ncRNA_CDS)

  #= RNAfold
  path_RNAfold = ut.find_RNAfold_path () #mirdup needs it
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
  nbTargets = paramDict['nbTargets']

  #= differentail analysis
  perform_differnatial_analysis = paramDict['perform_differnatial_analysis']
  diffguide_file = paramDict['diffguide_file']

  #= KEGG annotation
  gene_vs_pathway_file =  paramDict['gene_vs_pathway_file']

  #= enrichment analysis
  perform_KEGGpathways_enrichment_analysis= paramDict['perform_KEGGpathways_enrichment_analysis']
  pathway_description_file = paramDict['pathway_description_file']

  #= end of paramDict naming =================================================================================

  #= make required folders if not exist
  reps = [rep_output, rep_tmp]
  ut.makedirs_reps (reps)

  #= addFile
  sc.addPyFile(project_path + '/src/utils.py')
  sc.addPyFile(project_path + '/src/mirLibRules.py')
  sc.addFile(project_path + '/src/eval_mircheck.pl')
  sc.addFile(project_path + '/lib/miRcheck.pm')
  sc.addFile(project_path + '/lib/miRdup_1.4/lib/weka.jar')
  sc.addFile(mirdup_jar)
  sc.addFile(mirdup_model)
  sc.addFile(project_path + '/lib/VARNAv3-93.jar')

  #= Objects for rule functions
  dmask_obj = mru.prog_dustmasker()
  dmask_cmd, dmask_env = dmask_obj.dmask_pipe_cmd()
  kn_obj = mru.prog_knownNonMiRNA(broadcastVar_d_ncRNA_CDS.value)
  rnafold_obj = mru.prog_RNAfold(temperature)
  mircheck_obj = mru.prog_mirCheck(mcheck_param, project_path)
  mirdup_obj = mru.prog_miRdup (rep_tmp, mirdup_model, mirdup_jar, path_RNAfold)
  profile_obj = mru.prog_dominant_profile()
  miranda_obj = mru.prog_miRanda(Max_Score_cutoff, Max_Energy_cutoff, target_file, rep_tmp, miranda_binary, Gap_Penalty, nbTargets)


  #= Fetch library files in rep_input
  infiles = [f for f in listdir(rep_input) if os.path.isfile(os.path.join(rep_input, f))]
  print('============================================================\n')
  print('infiles:')
  for infile in infiles: print(infile)
  #= Time processing of libraries
  timeDict = {}

    
  print('\n====================== mirLibSpark =========================')
  print('====================== ' + appId + ' =================')
  print('============================================================\n')
  time_a = datetime.datetime.now()
  print(time_a, 'begin time')
  
  for infile in infiles :
    if infile[-1:] == '~': continue
    print ("--Processing of the library: ", infile)

    inBasename = os.path.splitext(infile)[0] #= lib name
    infile = rep_input+infile

    if input_type == 'q': #= fastq
      infile = ut.convert_fastq_file_to_KeyValue(infile, rep_tmp, inBasename)

    print ("  Start of miRNA prediction...", end="\n")
    print(datetime.datetime.now(), 'start')
    startLib = time.time()
    
    #= Convert the text file to RDD object
    ## in : file
    ## out: (a) u'seq\tfreq', 
    ##      (b) u'seq1', u'seq2', u'seq1', 
    ##      (c) u'>name1\nseq1', u'>name2\nseq2', u'>name3\nseq1',
    ##      (d) u'seq\tquality'
    distFile_rdd = sc.textFile("file:///" + infile, partition) #= partition is 2 if not set 
    if distFile_rdd.isEmpty():
      print(infile, 'is an empty file, omit this file')
      continue
    if reporting == 1: print(datetime.datetime.now(), 'NB distFile_rdd: ', distFile_rdd.count())#

    #= Unify different input formats to "seq freq" elements
    if input_type == 'w': #= raw
    ## in : u'seq\tfreq'
    ## out: ('seq', freq)
      ## note that type_a does not need to collapse nor trim.
      collapse_rdd = distFile_rdd.map(lambda line: mru.rearrange_rule(line, '\t'))
    else:
      if input_type == 'r': #= reads
      ## in : u'seq1', u'seq2', u'seq1'
      ## out: u'seq1', u'seq2', u'seq1'
        input_rdd = distFile_rdd

      elif input_type == 'a': #= fasta
      ## in : u'>name1\nseq1', u'>name2\nseq2', u'>name3\nseq1'
      ## out: u'seq1', u'seq2', u'seq1'
        input_rdd = distFile_rdd.filter(lambda line: not line[0] == '>')

      elif input_type == 'q': #= processed fastq
      ## in : u'seq\tquality'
      ## out: u'seq1', u'seq2', u'seq1'
        input_rdd = distFile_rdd.map(lambda word: word.split('\t')[0])

    #= trim adapters
      if not adapter == 'none':
        trim_adapter_rdd = input_rdd.map(lambda e: ut.trim_adapter (e, adapter))
      else: trim_adapter_rdd = input_rdd
      
    #= colapse seq and calculate frequency
      ## in : u'seq1', u'seq2', u'seq1'
      ## out: ('seq', freq)
      collapse_rdd = trim_adapter_rdd.map(lambda word: (word, 1))\
                                     .reduceByKey(lambda a, b: a+b)
  
    #= Filtering sRNA low frequency
    ## in : ('seq', freq)
    ## out: ('seq', freq)
    sr_low_rdd = collapse_rdd.filter(lambda e: int(e[1]) > limit_srna_freq)
    if reporting == 1: print(datetime.datetime.now(), 'NB sr_low_rdd: ', sr_low_rdd.count())

    
    #= Filtering short length
    ## in : ('seq', freq)
    ## out: ('seq', freq)
    sr_short_rdd = sr_low_rdd.filter(lambda e: len(e[0]) > limit_len).persist()  # TO KEEP IT, reused in bowFrq_rdd 
    if reporting == 1: print(datetime.datetime.now(), 'NB sr_short_rdd: ', sr_short_rdd.count())
    
    #= Filtering with DustMasker
    ## in : ('seq', freq)
    ## out: 'seq'
    dmask_rdd = sr_short_rdd.map(lambda e: '>s\n' + e[0])\
                            .pipe(dmask_cmd, dmask_env)\
                            .filter(lambda e: e.isupper() and not e.startswith('>'))\
                            .map(lambda e: str(e.rstrip()))\
                            .persist()
    if reporting == 1: print(datetime.datetime.now(), 'NB dmask_rdd: ', dmask_rdd.count())


    mergebowtie_rdd = sc.emptyRDD()
    for i in range(len(chromosomes)):
      ch = chromosomes[i]
      p = b_index_path + ch + '/' + bowtie_index_suffix
      bowtie_obj = mru.prog_bowtie(p)
      bowtie_cmd, bowtie_env = bowtie_obj.Bowtie_pipe_cmd()
      #================================================================================================================
      #================================================================================================================
      #= Mapping with Bowtie
      ## in : 'seq'
      ## out: ('seq', [nbLoc, [['strd','chr',posChr],..]])
      bowtie_rdd = dmask_rdd.map(lambda e: " "+e)\
                            .pipe(bowtie_cmd, bowtie_env)\
                            .map(bowtie_obj.bowtie_rearrange_map)\
                            .groupByKey()\
                            .map(lambda e: (e[0], [len(list(e[1])), list(e[1])]))
      mergebowtie_rdd = mergebowtie_rdd.union(bowtie_rdd)\
                                       .persist()
      #================================================================================================================
      #================================================================================================================
    if reporting == 1: print(datetime.datetime.now(), 'NB mergebowtie_rdd: ', mergebowtie_rdd.count())
    print(datetime.datetime.now(), 'mergebowtie_rdd')

    
    #= Getting the expression value for each reads
    ## in : ('seq', [nbLoc, [['strd','chr',posChr],..]])
    ## out: ('seq', [freq, nbLoc, [['strd','chr',posChr],..]])
    bowFrq_rdd = mergebowtie_rdd.join(sr_short_rdd)\
                                .map(bowtie_obj.bowtie_freq_rearrange_rule)\
                                .persist()
    if reporting == 1: print(datetime.datetime.now(), 'NB bowFrq_rdd: ', bowFrq_rdd.count())

    #= Filtering, keep miRNA length = 21, 22, 23, 24
    ## in : ('seq', [freq, nbLoc, [['strd','chr',posChr],..]])
    ## out: ('seq', [freq, nbLoc, [['strd','chr',posChr],..]])
    mr_meyers2018len_rdd = bowFrq_rdd.filter(lambda e: len(e[0]) < 25 and len(e[0]) > 20)
    if reporting == 1: print(datetime.datetime.now(), 'NB mr_meyers2018len_rdd: ', mr_meyers2018len_rdd.count())

    #= Filtering miRNA low frequency
    ## in : ('seq', [freq, nbLoc, [['strd','chr',posChr],..]])
    ## out: ('seq', [freq, nbLoc, [['strd','chr',posChr],..]])
    mr_low_rdd = mr_meyers2018len_rdd.filter(lambda e: e[1][0] > limit_mrna_freq)
    if reporting == 1: print(datetime.datetime.now(), 'NB mr_low_rdd: ', mr_low_rdd.count())
  
    #= Filtering high nbLocations and zero location
    ## in : ('seq', [freq, nbLoc, [['strd','chr',posChr],..]])
    ## out: ('seq', [freq, nbLoc, [['strd','chr',posChr],..]])
    nbLoc_rdd = mr_low_rdd.filter(lambda e: e[1][1] > 0 and e[1][1] < limit_nbLoc)
    if reporting == 1: print(datetime.datetime.now(), 'NB nbLoc_rdd: ', nbLoc_rdd.count())
    
    #= Flatmap the RDD
    ## in : ('seq', [freq, nbLoc, [['strd','chr',posChr],..]])
    ## out: ('seq', [freq, nbLoc, ['strd','chr',posChr])
    flat_rdd = nbLoc_rdd.flatMap(mru.flatmap_mappings)
    print(datetime.datetime.now(), 'flat_rdd')
    
    #= Filtering known non-miRNA ##
    ## in : ('seq', [freq, nbLoc, ['strd','chr',posChr])
    ## out: ('seq', [freq, nbLoc, ['strd','chr',posChr])
    excluKnownNon_rdd = flat_rdd.filter(kn_obj.knFilterByCoor)
    if reporting == 1: print(datetime.datetime.now(), 'excluKnownNon_rdd distinct: ', excluKnownNon_rdd.groupByKey().count())


    mergeChromosomesResults_rdd = sc.emptyRDD()
    for i in range(len(chromosomes)):
      ch = chromosomes[i]
      broadcastVar_genome = sc.broadcast(ut.getGenome(genome_path, ".fa", ch))
      prec_obj = mru.extract_precurosrs(broadcastVar_genome.value, pri_l_flank, pri_r_flank, pre_flank)
      #================================================================================================================
      #================================================================================================================
      #= Extraction of the pri-miRNA
      ## in : ('seq', [freq, nbLoc, ['strd','chr',posChr])
      ## out: ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri]])
      primir_rdd = excluKnownNon_rdd.filter(prec_obj.hasKey)\
                                    .flatMap(prec_obj.extract_prim_rule)
      mergeChromosomesResults_rdd = mergeChromosomesResults_rdd.union(primir_rdd)\
                                                               .persist()
      broadcastVar_genome.unpersist()
      #================================================================================================================
      #================================================================================================================
    if reporting == 1: print(datetime.datetime.now(), 'NB mergeChromosomesResults: ', mergeChromosomesResults_rdd.groupByKey().count())
    print(datetime.datetime.now(), 'mergeChromosomesResults_rdd') #= BOTTLE NECK
    
    #= pri-miRNA folding
    ## in : ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri]])
    ## out: ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold']])
    pri_fold_rdd = mergeChromosomesResults_rdd.map(lambda e: rnafold_obj.RNAfold_map_rule(e, 3))

    #= Validating pri-mirna with mircheck, keep unique precursors
    ## in : ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold']])
    ## out: ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold','mkPred','mkStart','mkStop']])
    pri_mircheck_rdd = pri_fold_rdd.map(lambda e: mircheck_obj.mirCheck_map_rule(e, 3))\
                                   .filter(lambda e: any(e[1][3]))\
                                   .map(lambda e: (e[0] + e[1][2][0] + e[1][2][1] + str(e[1][2][2]) + e[1][3][4] + e[1][3][5], e)  )\
                                   .reduceByKey(lambda a, b: a)\
                                   .map(lambda e: e[1])
    if reporting == 1: print(datetime.datetime.now(), 'NB pri_mircheck_rdd: ', pri_mircheck_rdd.groupByKey().count())
    print(datetime.datetime.now(), 'pri_mircheck_rdd') #= BOTTLE NECK


    #= Filtering len(pre-mirna) < 301 nt
    len300_rdd = pri_mircheck_rdd.filter(lambda e: (int(e[1][3][5]) - int(e[1][3][4])) < 301)
    if reporting == 1: print(datetime.datetime.now(), 'NB len300_rdd: ', len300_rdd.groupByKey().count())
  

    #======================#
    #= REPARTITION        =#
    #======================#
    #= Filtering structure with branched loop
    ## in : ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold','mkPred','mkStart','mkStop']])
    ## out: ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold','mkPred','mkStart','mkStop']])
    one_loop_rdd = len300_rdd.filter(lambda e: ut.containsOnlyOneLoop(e[1][3][2][int(e[1][3][4]) : int(e[1][3][5])+1]))\
                             .repartition(partition)
    if reporting == 1: print(datetime.datetime.now(), 'NB one_loop_rdd distinct : ', one_loop_rdd.groupByKey().count())



    #= Extraction of the pre-miRNA
    ## in : ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold','mkPred','mkStart','mkStop']])
    ## out: ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold','mkPred','mkStart','mkStop'], ['preSeq',posMirPre]])
    premir_rdd = one_loop_rdd.map(lambda e: prec_obj.extract_prem_rule(e, 3))

    
    #= pre-miRNA folding
    ## in : ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold', 'mkPred','mkStart','mkStop'], ['preSeq',posMirPre]])
    ## out: ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold', 'mkPred','mkStart','mkStop'], ['preSeq',posMirPre,'preFold']])
    pre_fold_rdd = premir_rdd.map(lambda e: rnafold_obj.RNAfold_map_rule(e, 4))
    #print('NB pre_fold_rdd: ', pre_fold_rdd.count())
    print(datetime.datetime.now(), 'pre_fold_rdd') 



    #================================================================================================================
    #= DEFUNCT: Validating pre-mirna with mircheck II -- replaced by mirdup
    ## in : ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold', 'mkPred','mkStart','mkStop'], ['preSeq',posMirPre,'preFold']])
    ## out: ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold', 'mkPred','mkStart','mkStop'], ['preSeq',posMirPre,'preFold','mkPred','mkStart','mkStop']])
    #pre_vld_rdd0 = pre_fold_rdd.map(lambda e: mircheck_obj.mirCheck_map_rule(e, 4))\
                              #.filter(lambda e: any(e[1][4]))
    #================================================================================================================
   

    #======================#
    #= REPARTITION        =#
    #======================#
    #= Validating pre-mirna with miRdup zipWithUniqueId
    ## in : ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold', 'mkPred','mkStart','mkStop'], ['preSeq',posMirPre,'preFold']])
    ## out: ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold', 'mkPred','mkStart','mkStop'], ['preSeq',posMirPre,'preFold','mpPred','mpScore']])
    pre_mirdup_rdd = pre_fold_rdd.zipWithIndex()\
                                 .map(mirdup_obj.run_miRdup)\
                                 .filter(lambda e: e[1][4][3] == "true")\
                                 .repartition(partition)
    if reporting == 1: print(datetime.datetime.now(), 'NB pre_mirdup_rdd distinct: ', pre_mirdup_rdd.groupByKey().count())
    print(datetime.datetime.now(), 'pre_mirdup_rdd distinct') #= BOTTLE NECK
    
    
    #= Filtering by expression profile (< 20%)
    ## in : ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold', 'mkPred','mkStart','mkStop'], ['preSeq',posMirPre,'preFold','mpPred','mpScore']])
    ## out: ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold', 'mkPred','mkStart','mkStop'], ['preSeq',posMirPre,'preFold','mpPred','mpScore'], totalfrq])

    #= keys = ['3A+', '5D-', ...]
    profile_keyvalue_rdd = pre_mirdup_rdd.map(lambda e: (e[1][2][1] + e[1][2][0], e))\
                                         .persist()
    keys_chromo_strand = profile_keyvalue_rdd.map(lambda e: (e[0], 1))\
                                             .reduceByKey(lambda a, b: a+b)\
                                             .map(lambda e: e[0])\
                                             .collect()
    mergeProfileChromo_rdd = sc.emptyRDD()
    x_rdd = bowFrq_rdd.flatMap(mru.flatmap_mappings)\
                      .map(lambda e: (e[1][2][1] + e[1][2][0], [e[1][2][2], e[1][0]]) )
    for chromo_strand in keys_chromo_strand:
      #================================================================================================================
      #================================================================================================================
      #======================#
      #= REPARTITION x2     =#
      #======================#
      y_rdd = x_rdd.filter(lambda e: e[0] == chromo_strand)
      broadcastVar_dict_bowtie_chromo_strand = sc.broadcast(profile_obj.get_bowtie_strandchromo_dict(y_rdd.collect()))
      profile_value_rdd = profile_keyvalue_rdd.filter(lambda e: e[0] == chromo_strand)\
                                              .repartition(partition)\
                                              .map(lambda e: profile_obj.computeProfileFrq(e[1], broadcastVar_dict_bowtie_chromo_strand.value))\
                                              .filter(lambda e: e[1][0] / (float(e[1][5]) + 0.1) > 0.2)
      mergeProfileChromo_rdd = mergeProfileChromo_rdd.union(profile_value_rdd)\
                                                     .repartition(partition)\
                                                     .persist()
      #================================================================================================================
      #================================================================================================================

    slim_rdd = mergeProfileChromo_rdd.map(mru.slimrule)


    if reporting == 1: print(datetime.datetime.now(), 'NB slim_rdd NON distinct: ', slim_rdd.count())
    print('NB slim_rdd distinct: ', slim_rdd.groupByKey().count()) #= always report the nb of final prediction
    print(datetime.datetime.now(), 'slim_rdd')
    
    #= collecting final miRNA predictions
    libresults = slim_rdd.collect()
    print(datetime.datetime.now(), 'libresults=slim_rdd.collect()')#= BOTTLE NECK

    endLib = time.time() 
    timeDict[inBasename] = endLib - startLib
    print ('miRNA predcition time for lib ', inBasename, ': ', timeDict[inBasename])
    print ('  End of miRNA prediction     ', end='\n')

    #= write results to a file
    eachLiboutFile = rep_output  +  appId + '_miRNAprediction_' + inBasename + '.txt'
    ut.writeToFile (libresults, eachLiboutFile)

  #= print executions time  to a file
  outTime = rep_output + appId + '_time.txt'
  ut.writeTimeLibToFile (timeDict, outTime, appId, paramDict)

  
  #===============================================================================================================
  #=
  #= Post processing after miRNA prediction
  #=
  #===============================================================================================================
  #= make summary table of all libraries in one submission with expressions in the field
  keyword = appId + '_miRNAprediction_'
  infiles = [f for f in listdir(rep_output) if (os.path.isfile(os.path.join(rep_output, f)) and f.startswith(keyword))]
  Precursor, distResultSmallRNA = ut.writeSummaryExpressionToFile (infiles, rep_output, appId)

  ### in : ( 'seq' )
  ### out: ( 'seq', zipindex)
  distResultSmallRNA_rdd = sc.parallelize(distResultSmallRNA, partition)\
                             .zipWithIndex()
  

  #= varna precursor visualization
  varna_obj = mru.prog_varna(appId, rep_output) 
  ## in : ([miRNAseq, strand, chromo, posChr, preSeq, posMirPre, preFold, mkPred, newfbstart, newfbstop, mpPred, mpScore], zipindex)
  Precursor_rdd = sc.parallelize(Precursor, partition)\
                    .zipWithIndex()
  distResultSmallRNA = distResultSmallRNA_rdd.collect()
  d_rna_index = {} # {seq: index}
  for i in distResultSmallRNA: d_rna_index[ i[0] ] = i[1]
  ## out : ( PrecursorIndex, miRNAseq, strand, chromo, posChr, preSeq, posMirPre, preFold, mkPred, newfbstart, newfbstop, mpPred, mpScore, miRNAindex )
  PrecursorVis = Precursor_rdd.map(varna_obj.run_VARNA)\
                              .map(lambda e: mru.matchRNAidRule(e, d_rna_index))\
                              .collect()
  ut.write_index (PrecursorVis, rep_output, appId)
  print('PrecursorVis done')
  
  
  #= miranda target prediction
  ## in : ('miRNAseq', zipindex)
  ## out: ('miRNAseq', [[target1 and its scores], [target2 and its scores]])
  sc.addFile(miranda_binary)
  sc.addFile(target_file)
  mirna_and_targets = distResultSmallRNA_rdd.map(miranda_obj.computeTargetbyMiranda)\
                                            .collect()
  ut.writeTargetsToFile (mirna_and_targets, rep_output, appId)
  print('Target prediction done')
  

  ## I dont know what is the use of this, maybe there is no use...
  ## in: ('miRNAseq', [[targetgene1 and its scores], [targetgene2 and its scores]])
  ## out:( 'targetgene' )
  #master_distinctTG = miranda_rdd.map(lambda e: [  i[0].split('.')[0] for i in e[1]  ])\
  #                               .reduce(lambda a, b: a+b)
  #master_distinctTG = sorted(list(set(master_distinctTG)))
  #print( master_distinctTG )


  #= clear caches (memory leak)
  broadcastVar_paramDict.unpersist()
  dmask_rdd.unpersist()
  sr_short_rdd.unpersist()
  mergebowtie_rdd.unpersist()
  mergeChromosomesResults_rdd.unpersist()
  broadcastVar_d_ncRNA_CDS.unpersist()
  bowFrq_rdd.unpersist()
  mergeProfileChromo_rdd.unpersist()
  broadcastVar_dict_bowtie_chromo_strand.unpersist()
  profile_keyvalue_rdd.unpersist()


  #= end of spark context, stop to allow running multiple SparkContexts
  sc.stop() 
  print(datetime.datetime.now(), 'sc stop time')
  #===============================================================================================================
  #===============================================================================================================
  ##appId = 'local-1538502520294'
  #= diff analysis 
  if perform_differnatial_analysis == 'yes':
    diff_outs = ut.diff_output(diffguide_file, rep_output, appId)
    print('Differential analysis done')

  if perform_KEGGpathways_enrichment_analysis == 'yes':
    #= KEGG annotation
    list_mirna_and_topscoredTargetsKEGGpathway = ut.annotate_target_genes_with_KEGGpathway (gene_vs_pathway_file, rep_output, appId)
    print('KEGG pathway annotation done')
    #= KEGG enrichment analysis 
    ut.perform_enrichment_analysis (diff_outs, pathway_description_file, list_mirna_and_topscoredTargetsKEGGpathway, rep_output, appId, project_path)
    print('\nKEGG pathway enrichment analysis done')
  #===============================================================================================================
  #===============================================================================================================


  time_b = datetime.datetime.now()
  print(time_b, 'finish time')
  print('total running time: ', time_b - time_a)
  print('====================== End of ' + appId + ' =============\n')






