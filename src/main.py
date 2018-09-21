import os

#cmd = 'python purge.py'
#os.system (cmd)

cmd = 'git pull origin join'
os.system (cmd)



#cmd = 'cp ../input_samples/fake_a3.txt ../input'
#cmd = 'rm -fr ../input/fake_a3.txt'
#os.system(cmd)


#cmd = 'cp ../input_samples/100.txt ../input'
cmd = 'rm -fr ../input/100.txt'
os.system(cmd)





cmd = 'time spark-submit mirLibPipeline.py ../paramfile.txt 2>/dev/null'
os.system(cmd)


#cmd = 'mkdir ../dbs/ATH/bowtie_split_index'
#os.system(cmd)


#chro = [1, 2, 3, 4, 5, 'C', 'M']
#for c in chro:
#  refgenome = '../dbs/ATH/Genome/TAIR10_chr' + str(c) + '.fas'
#  cmd = '../lib/bowtie-build -f ' + refgenome + ' atht10_' + 'chr' + str(c)
#  os.system(cmd)



#chro = [1, 2, 3, 4, 5, 'C', 'M']
#path = '../dbs/ATH/bowtie_split_index/chr'
#key = '/atht10_chr'
#newkey = '/a_thaliana_t10_chr'
#for c in chro:
#  prefix = path + str(c) + key + str(c)
#  newprefix = path + str(c) + newkey + str(c)
#  s1 = '.1.ebwt'
#  s2 = '.2.ebwt'
#  s3 = '.3.ebwt'
#  s4 = '.4.ebwt'
#  s5 = '.rev.1.ebwt'
#  s6 = '.rev.2.ebwt'
#  suffixes = [s1, s2, s3, s4, s5, s6]
#  for s in suffixes:
#    cmd = 'git mv ' + prefix + s + ' ' + newprefix + s
#    os.system(cmd)



#chro = [0]
#path = '../dbs/ATH/bowtie_split_index/chr'
#key = '/a_thaliana_t10'
#for c in chro:
#  prefix = path + str(c) + key
#  newprefix = path + 'All' + key
#  s1 = '.1.ebwt'
#  s2 = '.2.ebwt'
#  s3 = '.3.ebwt'
#  s4 = '.4.ebwt'
#  s5 = '.rev.1.ebwt'
#  s6 = '.rev.2.ebwt'
#  suffixes = [s1, s2, s3, s4, s5, s6]
#  for s in suffixes:
#    cmd = 'git mv ' + prefix + s + ' ' + newprefix + s
#    os.system(cmd)




