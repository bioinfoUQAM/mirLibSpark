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


if __name__ == '__main__' :

  if not len(sys.argv) == 2:
    sys.stderr.write('Three arguments required\nUsage: spark-submit mirLibPipeline.py <path to paramfile> 2>/dev/null\n')
    sys.exit()


  paramfile = sys.argv[1]
  paramDict = ut.readParam (paramfile)
  #= EXMAMINE OPTIONS 
  print('\nVerifying parameters ...')
  ut.validate_options(paramDict)


  #= spark configuration
  appMaster = paramDict['sc_master']                #"local[*]" 
  appName = paramDict['sc_appname']                 #"mirLibSpark"
  mstrMemory = paramDict['sc_mstrmemory']           #"4g"
  execMemory = paramDict['sc_execmemory']           #"4g"
  execCores = paramDict['sc_execcores']             #2
  partition = int(paramDict['sc_partition'])

  #= Spark context
  sc = ut.pyspark_configuration(appMaster, appName, mstrMemory, execMemory, execCores)

  #= Spark application ID
  appId = str(sc.applicationId)

  #= broadcast paramDict
  broadcastVar_paramDict = sc.broadcast(paramDict)
  paramDict = broadcastVar_paramDict.value

  #= Parameters and cutoffs =========================
  #= paths
  input_type = paramDict['input_type']
  adapter = ut.tr_U_T (paramDict['adapter'])
  project_path = paramDict['project_path'][:-1]
  rep_input = paramDict['input_path']
  rep_output = paramDict['output_path']
  rep_msub_jobsOut = project_path + '/workdir/jobsOut'
  rep_tmp = project_path + '/tmp/'     

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
  reps = [rep_output, rep_tmp, rep_msub_jobsOut]
  ut.makedirs_reps (reps)

  #= addFile
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

  #= Objects for rule functions
  dmask_obj = mru.prog_dustmasker()
  dmask_cmd, dmask_env = dmask_obj.dmask_pipe_cmd()
  kn_obj = mru.prog_knownNonMiRNA(broadcastVar_d_ncRNA_CDS.value)
  rnafold_obj = mru.prog_RNAfold(temperature)
  mircheck_obj = mru.prog_mirCheck(mcheck_param)
  mirdup_obj = mru.prog_miRdup (rep_tmp, mirdup_model, mirdup_jar, path_RNAfold)
  profile_obj = mru.prog_dominant_profile()
  miranda_obj = mru.prog_miRanda(Max_Score_cutoff, Max_Energy_cutoff, target_file, rep_tmp, miranda_binary, Gap_Penalty, nbTargets)


  #= Fetch library files in rep_input
  infiles = [f for f in listdir(rep_input) if os.path.isfile(os.path.join(rep_input, f))]

  #= Time processing of libraries
  timeDict = {}

  #== test checkpoint()
  #sc.setCheckpointDir(rep_output)
    
  print('\n====================== mirLibSpark =========================')
  print('====================== ' + appId + ' =================')
  #for k, v in paramDict.items(): print(k, ': ', v)
  print('============================================================\n')
  print('begin time:', datetime.datetime.now())
  #'''
  libRESULTS = [] 
  for infile in infiles :
    if infile[-1:] == '~': continue
    print ("--Processing of the library: ", infile)

    inBasename = os.path.splitext(infile)[0] #= lib name
    infile = rep_input+infile
    inKvfile = rep_tmp + inBasename + '.kv.txt'
    # hdfsFile = inBasename + '.hkv.txt' #= spark will take care of textFile() and pass it through hdfs, so no need this line

    if input_type == 'd': #= fastq
      ut.convert_fastq_file_to_KeyValue(infile, inKvfile)
      infile = inKvfile
      
    print ("  Start of the processing...", end="\n")
    startLib = time.time()
    
    #= Convert the text file to RDD object
    ## in : file
    ## out: (a) u'seq\tfreq', 
    ##      (b) u'seq1', u'seq2', u'seq1', 
    ##      (c) u'>name1\nseq1', u'>name2\nseq2', u'>name3\nseq1',
    ##      (d) u'seq\tquality'
    distFile_rdd = sc.textFile("file:///" + infile, partition) #= partition is 2 if not set 
    #print('NB distFile_rdd: ', distFile_rdd.count())#

    #= Unify different input formats to "seq freq" elements
    if input_type == 'a': #= raw
    ## in : u'seq\tfreq'
    ## out: ('seq', freq)
      ## note that type_a does not need to collapse nor trim.
      collapse_rdd = distFile_rdd.map(lambda line: mru.rearrange_rule(line, '\t'))#\
    else:
      if input_type == 'b': #= reads
      ## in : u'seq1', u'seq2', u'seq1'
      ## out: u'seq1', u'seq2', u'seq1'
        input_rdd = distFile_rdd

      elif input_type == 'c': #= fasta
      ## in : u'>name1\nseq1', u'>name2\nseq2', u'>name3\nseq1'
      ## out: u'seq1', u'seq2', u'seq1'
        input_rdd = distFile_rdd.filter(lambda line: not line[0] == '>')

      elif input_type == 'd': #= processed fastq
      ## in : u'seq\tquality'
      ## out: u'seq1', u'seq2', u'seq1'
        input_rdd = distFile_rdd.map(lambda word: word.split('\t')[0])

    #= trim adapters
      if not adapter == 'none':
        trim_adapter_rdd = input_rdd.map(lambda e: trim_adapter (e, adapter))
      else: trim_adapter_rdd = input_rdd
      
    #= colapse seq and calculate frequency
      ## in : u'seq1', u'seq2', u'seq1'
      ## out: ('seq', freq)
      collapse_rdd = trim_adapter_rdd.map(lambda word: (word, 1)).reduceByKey(lambda a, b: a+b)
  
    #= Filtering sRNA low frequency
    ## in : ('seq', freq)
    ## out: ('seq', freq)
    sr_low_rdd = collapse_rdd.filter(lambda e: int(e[1]) > limit_srna_freq)
    #print('NB sr_low_rdd: ', sr_low_rdd.count())####################################################
    
    #= Filtering short length
    ## in : ('seq', freq)
    ## out: ('seq', freq)
    sr_short_rdd = sr_low_rdd.filter(lambda e: len(e[0]) > limit_len).persist()  # TO KEEP IT
    #print('NB sr_short_rdd: ', sr_short_rdd.count())###################################################
    
    #= Filtering with DustMasker
    ## in : ('seq', freq)
    ## out: 'seq'
    dmask_rdd = sr_short_rdd.map(lambda e: '>s\n'+e[0])\
                            .pipe(dmask_cmd, dmask_env)\
                            .filter(lambda e: e.isupper() and not e.startswith('>'))\
                            .map(lambda e: str(e.rstrip()))\
                            .persist()
    #print('NB dmask_rdd: ', dmask_rdd.count())############################################

    mergebowtie_rdd = sc.emptyRDD()
    for i in range(len(chromosomes)):
      ch = chromosomes[i]
      #= furture work: case insensitive
      p = b_index_path + ch.replace('Chr', 'chr') + '/' + bowtie_index_suffix + '_' + ch.replace('chr', 'Chr') 
      if ch == 'All': p = b_index_path + ch + '/' + bowtie_index_suffix 
      bowtie_obj = mru.prog_bowtie(p)
      bowtie_cmd, bowtie_env = bowtie_obj.Bowtie_pipe_cmd()
      #================================================================================================================
      #================================================================================================================
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
      #print('NB bowtie_rdd: ', len(bowtie_rdd.collect()))##################################################
      #================================================================================================================
      #================================================================================================================
      #================================================================================================================
      #================================================================================================================
      mergebowtie_rdd = mergebowtie_rdd.union(bowtie_rdd).persist()

    #= Getting the expression value for each reads
    ## in : ('seq', [nbLoc, [['strd','chr',posChr],..]])
    ## out: ('seq', [freq, nbLoc, [['strd','chr',posChr],..]])
    bowFrq_rdd = mergebowtie_rdd.join(sr_short_rdd)\
                           .map(bowtie_obj.bowtie_freq_rearrange_rule)
    #print('NB bowFrq_rdd: ', bowFrq_rdd.count())############################################ 
    #180921 fake_a.txt takes 17 secs till this step, option chromo=All
    #180921 100.txt takes 23 secs till this step, option chromo=All

    #= Create dict, chromo_strand as key to search bowtie blocs in the following dict
    dict_bowtie_chromo_strand = profile_obj.get_bowtie_strandchromo_dict(bowFrq_rdd.collect())
    broadcastVar_bowtie_chromo_strand = sc.broadcast(dict_bowtie_chromo_strand) #= get the value by broadcastVar.value


    #= Filtering miRNA low frequency
    ## in : ('seq', [freq, nbLoc, [['strd','chr',posChr],..]])
    ## out: ('seq', [freq, nbLoc, [['strd','chr',posChr],..]])
    mr_low_rdd = bowFrq_rdd.filter(lambda e: e[1][0] > limit_mrna_freq)
    #print('NB mr_low_rdd: ', mr_low_rdd.count())##################################################
    
    #= Filtering high nbLocations and zero location
    ## in : ('seq', [freq, nbLoc, [['strd','chr',posChr],..]])
    ## out: ('seq', [freq, nbLoc, [['strd','chr',posChr],..]])
    nbLoc_rdd = mr_low_rdd.filter(lambda e: e[1][1] > 0 and e[1][1] < limit_nbLoc)
    #print('NB nbLoc_rdd: ', nbLoc_rdd.count())###################################################
    
    
    #= Flatmap the RDD
    ## in : ('seq', [freq, nbLoc, [['strd','chr',posChr],..]])
    ## out: ('seq', [freq, nbLoc, ['strd','chr',posChr])
    flat_rdd = nbLoc_rdd.flatMap(mru.flatmap_mappings)
    #print('NB flat_rdd distinct (this step flats elements): ', flat_rdd.groupByKey().count())##################
    #print('NB flat_rdd not distinct: ', flat_rdd.count())##################

    ###############################
    ## Filtering known non-miRNA ##
    ###############################
    ## in : ('seq', [freq, nbLoc, ['strd','chr',posChr])
    ## out: ('seq', [freq, nbLoc, ['strd','chr',posChr])
    #excluKnownNon_rdd = flat_rdd.repartition(100).filter(kn_obj.knFilterByCoor)
    excluKnownNon_rdd = flat_rdd.filter(kn_obj.knFilterByCoor)
    #print('excluKnownNon_rdd distinct: ', excluKnownNon_rdd.groupByKey().count())########
    

    mergeChromosomesResults_rdd = sc.emptyRDD()
    for i in range(len(chromosomes)):
      ch = chromosomes[i].replace('chr', 'Chr')
      genome = ut.getGenome(genome_path, ".fas", ch)
      broadcastVar_genome = sc.broadcast(genome)
      v = broadcastVar_genome.value
      prec_obj = mru.extract_precurosrs(v, pri_l_flank, pri_r_flank, pre_flank, ch)
      #================================================================================================================
      #================================================================================================================
      #================================================================================================================
      #================================================================================================================
      #= Extraction of the pri-miRNA
      ## in : ('seq', [freq, nbLoc, ['strd','chr',posChr])
      ## out: ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri]])
      primir_rdd = excluKnownNon_rdd.filter(prec_obj.hasKey)\
                                    .flatMap(prec_obj.extract_prim_rule)

      #= pri-miRNA folding
      ## in : ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri]])
      ## out: ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold']])
      pri_fold_rdd = primir_rdd.map(lambda e: rnafold_obj.RNAfold_map_rule(e, 3))
    
      #= Validating pri-mirna with mircheck
      ## in : ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold']])
      ## out: ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold','mkPred','mkStart','mkStop']])
      pri_vld_rdd = pri_fold_rdd.map(lambda e: mircheck_obj.mirCheck_map_rule(e, 3))\
                                .filter(lambda e: any(e[1][3]))
      #print('NB pri_vld_rdd distinct (mircheck): ', pri_vld_rdd.groupByKey().count())##

      #= Filtering structure with branched loop
      ## in : ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold','mkPred','mkStart','mkStop']])
      ## out: ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold','mkPred','mkStart','mkStop']])
      one_loop_rdd = pri_vld_rdd.filter(lambda e: ut.containsOnlyOneLoop(e[1][3][2][int(e[1][3][4]) : int(e[1][3][5])+1]))
      #print('NB one_loop_rdd distinct : ', len(one_loop_rdd.groupByKey().collect()))##
    
      #= Extraction of the pre-miRNA
      ## in : ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold','mkPred','mkStart','mkStop']])
      ## out: ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold','mkPred','mkStart','mkStop'], ['preSeq',posMirPre]])
      premir_rdd = one_loop_rdd.map(lambda e: prec_obj.extract_prem_rule(e, 3)) ## use one-loop rule
      #premir_rdd = pri_vld_rdd.map(lambda e: prec_obj.extract_prem_rule(e, 3)) ## ignore one-loop rule
      #================================================================================================================
      #================================================================================================================
      #================================================================================================================
      #================================================================================================================
      mergeChromosomesResults_rdd = mergeChromosomesResults_rdd.union(premir_rdd).persist()#.checkpoint()
      broadcastVar_genome.unpersist()
    #print('mergeChromosomesResults: ', mergeChromosomesResults_rdd.count())######## 
    #180921 fake_a.txt takes 42 secs to run till this line (All chromo)
    #180921 fake_a.txt takes 307 secs to run till this line (split chromo)

    
    #= pre-miRNA folding
    ## in : ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold', 'mkPred','mkStart','mkStop'], ['preSeq',posMirPre]])
    ## out: ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold', 'mkPred','mkStart','mkStop'], ['preSeq',posMirPre,'preFold']])
    pre_fold_rdd = mergeChromosomesResults_rdd.map(lambda e: rnafold_obj.RNAfold_map_rule(e, 4))

    #= Validating pre-mirna with mircheck II -- replaced by mirdup
    ## in : ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold', 'mkPred','mkStart','mkStop'], ['preSeq',posMirPre,'preFold']])
    ## out: ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold', 'mkPred','mkStart','mkStop'], ['preSeq',posMirPre,'preFold','mkPred','mkStart','mkStop']])
    #pre_vld_rdd0 = pre_fold_rdd.map(lambda e: mircheck_obj.mirCheck_map_rule(e, 4))\
                              #.filter(lambda e: any(e[1][4]))
    #print('NB pre_vld_rdd0 distinct (mircheck II): (', pre_vld_rdd0.groupByKey().count()), ')')

    #= Validating pre-mirna with miRdup zipWithUniqueId
    ## in : ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold', 'mkPred','mkStart','mkStop'], ['preSeq',posMirPre,'preFold']])
    ## out: ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold', 'mkPred','mkStart','mkStop'], ['preSeq',posMirPre,'preFold','mpPred','mpScore']])
    pre_vld_rdd = pre_fold_rdd.zipWithIndex()\
                              .map(mirdup_obj.run_miRdup)\
                              .filter(lambda e: e[1][4][3] == "true")
    #print('NB pre_vld_rdd distinct (mirdup): ', pre_vld_rdd.groupByKey().count())##

    
    #= Filtering by expression profile (< 20%)
    ## in : ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold', 'mkPred','mkStart','mkStop'], ['preSeq',posMirPre,'preFold','mpPred','mpScore']])
    ## out: ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold', 'mkPred','mkStart','mkStop'], ['preSeq',posMirPre,'preFold','mpPred','mpScore'], totalfrq])
    profile_rdd = pre_vld_rdd.map(lambda e: profile_obj.computeProfileFrq(e, broadcastVar_bowtie_chromo_strand.value))\
                             .filter(lambda e: e[1][0] / float(e[1][5]) > 0.2)
    print('NB profile_rdd distinct: ', profile_rdd.groupByKey().count())##
    libresults = profile_rdd.collect()
    libRESULTS.append( [inBasename, libresults] )
   
    endLib = time.time() 
    timeDict[inBasename] = endLib - startLib
    print ("  End of the processing     ", end="\n")

    #= write results to a file
    eachLiboutFile = rep_output  +  appId + '_miRNAprediction_' + inBasename + '.txt'
    ut.writeToFile (libresults, eachLiboutFile)

  #= print executions time  to a file
  outTime = rep_output + appId + '_time.txt'
  ut.writeTimeLibToFile (timeDict, outTime, appId, paramDict)

  #= make summary table of all libraries in one submission with expressions in the field
  keyword = appId + '_miRNAprediction_'
  infiles = [f for f in listdir(rep_output) if (os.path.isfile(os.path.join(rep_output, f)) and f.startswith(keyword))]
  Precursor = ut.writeSummaryExpressionToFile (infiles, rep_output, appId)
  #broadcastVar_Precursor = sc.broadcast(Precursor)   

  ## in:  ( lib, ('seq', [...]) )
  ## out: ( 'seq', [...] )
  #broadcastVar_libRESULTS = sc.broadcast(libRESULTS)  
  #libRESULTS_rdd = sc.parallelize(broadcastVar_libRESULTS.value, partition).flatMap(lambda e: e[1]) 
  libRESULTS_rdd = sc.parallelize(libRESULTS, partition).flatMap(lambda e: e[1]) 

  ## in:  ( 'seq', [...] )
  ## out: ( 'seq' ) 
  master_predicted_distinctMiRNAs_rdd = libRESULTS_rdd.map(lambda e: e[0]).distinct()
  
  ## in:  ( 'seq' ) 
  ## out: ('miRNAseq', zipindex)
  distResultSmallRNA_rdd = master_predicted_distinctMiRNAs_rdd.zipWithIndex() 

  
  #= varna
  varna_obj = mru.prog_varna(appId, rep_output) 

  ## in:  ( 'seq', [...] ) 
  ## mid: [miRNAseq, frq, nbLoc, strand, chromo, posChr, mkPred, mkStart, mkStop, preSeq, posMirPre, newfbstart, newfbstop, preFold, mpPred, mpScore, totalFrq] 
  ## out : [miRNAseq, strand, chromo, posChr, preSeq, posMirPre, preFold, mkPred, newfbstart, newfbstop, mpPred, mpScore]
  #Precursor_rdd = sc.parallelize(broadcastVar_Precursor.value, partition)
  Precursor_rdd = sc.parallelize(Precursor, partition)
  
  ## in : [miRNAseq, strand, chromo, posChr, preSeq, posMirPre, preFold, mkPred, newfbstart, newfbstop, mpPred, mpScore]
  ## out : ([miRNAseq, strand, chromo, posChr, preSeq, posMirPre, preFold, mkPred, newfbstart, newfbstop, mpPred, mpScore], zipindex)
  VARNA_rdd = Precursor_rdd.zipWithIndex()\
                               .map(varna_obj.run_VARNA)
  indexVis = VARNA_rdd.collect()
  ut.write_index (indexVis, rep_output, appId)
  
  
  #= miranda
  ## in : ('miRNAseq', zipindex)
  ## out: ('miRNAseq', [[target1 and its scores], [target2 and its scores]])
  miranda_rdd = distResultSmallRNA_rdd.map(miranda_obj.computeTargetbyMiranda)
  mirna_and_targets = miranda_rdd.collect()
  ut.writeTargetsToFile (mirna_and_targets, rep_output, appId)
  #'''

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
  broadcastVar_bowtie_chromo_strand.unpersist()

  #'''

  #= end of spark context
  #broadcastVar_paramDict.unpersist()
  #broadcastVar_d_ncRNA_CDS.unpersist()

  sc.stop() #= allow to run multiple SparkContexts

  #===============================================================================================================
  #===============================================================================================================
  #===============================================================================================================
  #===============================================================================================================
  #appId = 'local-1538110614002'

  #= diff analysis 
  if perform_differnatial_analysis == 'yes':
    diffguide, _ = ut.read_diffguide(diffguide_file)
    diff_outs = ut.diff_output(diffguide, rep_output, appId)


  if perform_KEGGpathways_enrichment_analysis == 'yes':
    #= KEGG annotation
    list_mirna_and_topscoredTargetsKEGGpathway = ut.annotate_target_genes_with_KEGGpathway (gene_vs_pathway_file, rep_output, appId)

    #= KEGG enrichment analysis 
    ut.input_for_enrichment_analysis (diff_outs, pathway_description_file, list_mirna_and_topscoredTargetsKEGGpathway, rep_output, appId)


  #===============================================================================================================
  #===============================================================================================================
  #===============================================================================================================
  #===============================================================================================================



  print('finish time:', datetime.datetime.now())
  print('====================== End of ' + appId + ' =============\n')


  #os.system('rm -fr ' + rep_tmp)



