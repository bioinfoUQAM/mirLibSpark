'''
program: utils.py
author: Chao-Jung Wu
author: M.A.Remita
date: 2017-03-25
version: 0.00.01
'''


# Configure a spark context
def pyspark_configuration(appMaster, appName, appMemory):
  from pyspark import SparkConf, SparkContext
  conf = SparkConf()
  conf.setMaster(appMaster)
  conf.setAppName(appName)
  conf.set("spark.executor.memory", appMemory)
  return SparkContext(conf = conf)

# Convert a file to hadoop file
def convertTOhadoop(rfile, hdfsFile):
  import os
  print('if pre-existing in hdfs, the file would be deleted before the re-distribution of a new file with the same name.\n')
  os.system('hadoop fs -rm ' + hdfsFile) # force delete any pre-existing file in hdfs with the same name.
  os.system('hadoop fs -copyFromLocal ' + rfile + ' ' + hdfsFile)

# Convert a fasta file into a key value file
def covert_fasta_to_KeyValue(infile, outfile):
  fh = open (infile, 'r')
  DATA = fh.readlines()
  fh.close()
  dict_sRNA = {}
  for i in xrange(0, len(DATA), 2):
      key = DATA[i].rstrip('\n')[1:]
      value = DATA[i+1].rstrip('\n')
      dict_sRNA[key] = value
  fh_out = open (outfile, 'w')
  for k, v in dict_sRNA.items():
      print >>fh_out, k + ':' + v
  fh_out.close()


# Convert a seq abundance file into a key value file
def convert_seq_freq_file_to_KeyValue(infile, outfile, v_sep):
  fh = open (infile, 'r')
  fh_out = open (outfile, 'w')

  dict_sRNA = {}
  i = 1
  
  for line in fh:
    data = line.rstrip('\n').split('\t')
    value = data[0] + v_sep + data[1]
    
    # check if the read was treated before (redundancy)
    if data[0] not in dict_sRNA:
      key = str(i).zfill(9)
      i+= 1
      dict_sRNA[data[0]] = 1
      print >>fh_out, key + v_sep + value
      
  fh.close()
  fh_out.close()
  
  def getRevComp (seq):
    from string import maketrans
    
    intab = "ACGT"
    outab = "TGCA"
    trantab = maketrans(intab, outab)
    n_seq = seq.translate(trantab)
    return n_seq[::-1];
  
  def tr_T_U (seq):
    from string import maketrans
    trantab = maketrans("T", "U")
    return seq.translate(trantab)
  
  def tr_U_T (seq):
    from string import maketrans
    trantab = maketrans("U", "T")
    return seq.translate(trantab)
