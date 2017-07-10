'''
script for test solely mirdup in spark
first test in Cloudera

Author: Julie Wu
Date: 2017-07-10

'''
from pyspark import SparkContext 
import subprocess as sbp

def mirdup_rule (e):
    '''
    java -jar miRdup_1.4/miRdup.jar -v miRdup_1.4/testFiles/julie_sequencesToValidate2.txt -c miRdup_1.4//model/Viridiplantae.model -r /usr/local/bin/
    '''
    tmp_file = '/home/cloudera/workspace/mirLibHadoop/test/tmp_mirna_seq_test_mirdup.txt'
    e = e.split('\t')
    with open (tmp_file, 'w') as fh_tmp:
      print >> fh_tmp, 'seqx\t' + e[1] + '\t' + e[2] #name, mirna, pre-mirna
    
    cmd = ['java', '-jar', '../lib/miRdup_1.4/miRdup.jar', '-v', tmp_file, '-c', '../lib/miRdup_1.4//model/Viridiplantae.model', '-r', '/usr/local/bin/']
    sproc = sbp.Popen(cmd, stdout=sbp.PIPE, shell=False)
    mirdupout = sproc.communicate()[0].split('\n')

    for i in mirdupout:
      data = i.split()
      if len(data) == 6 and data[0] == 'Correctly':
        mirdup_verdict = data[3]
        if mirdup_verdict == '1':
          return True
    return False

    if mirdup_verdict == '0':
      return False
    return True

def test_mirdup ():
    '''
    cloudera [~/workspace/miRNA_predictor/src]
    $ spark-submit 170710_test_mirdup_in_spark.py 2>/dev/null
    '''
    infile = 'file:///home/cloudera/workspace/mirLibHadoop/test/julie_sequencesToValidate_2T2F.txt'
    sc = SparkContext("local", "App Name")
    rdd1 = sc.textFile(infile)
    #rdd2 = rdd1.map(mirdup_rule)
    rdd2 = rdd1.filter(mirdup_rule)
    print(rdd2.collect())
    print('test_mirdup () ends')
 
test_mirdup ()
