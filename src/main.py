import os



#cmd = 'python purge.py';os.system (cmd)

cmd = 'git pull origin testWheatInHouse';os.system (cmd)

cmd = 'git branch -D ensembl';os.system(cmd)
cmd = 'git remote prune origin';os.system(cmd)
#cmd = 'git branch';os.system(cmd)
#cmd = 'git fetch origin';os.system(cmd)
#cmd = 'git checkout master';os.system(cmd)


#cmd = 'cp ../input_samples/fake_a.txt ../input';os.system(cmd)
cmd = 'rm -fr ../input/fake_a.txt';os.system(cmd)

#cmd = 'cp ../input_samples/fake_a3.txt ../input';os.system(cmd)
cmd = 'rm -fr ../input/fake_a3.txt';os.system(cmd)

#cmd = 'cp ../input_samples/fake_a5.txt ../input';os.system(cmd)
cmd = 'rm -fr ../input/fake_a5.txt';os.system(cmd)

#cmd = 'cp ../input_samples/100.txt ../input';os.system(cmd)
cmd = 'rm -f ../input/100.txt';os.system(cmd)

#cmd = 'mv ../input_samples/11w2013_t2_1.fasta ../input';os.system(cmd)
#cmd = 'mv ../input/11w2013_t2_1.fasta ../input_samples';os.system(cmd)



cmd = 'cp ../input_samples/wheat_a1.txt.txt ../input';os.system(cmd)
#cmd = 'rm -fr ../input/wheat_a1.txt';os.system(cmd)

cmd = 'cp ../input_samples/wheat_a4.txt.txt ../input';os.system(cmd)
#cmd = 'rm -fr ../input/wheat_a4.txt';os.system(cmd)



#cmd = 'git add -f ../output/*';os.system(cmd)
#cmd = 'git commit -m "save a copy of output"';os.system(cmd)
#cmd = 'git push origin refactor';os.system(cmd)



cmd = 'free -m';os.system(cmd)

#########################################################################################
#cmd = 'time spark-submit mirLibPipeline.py ../paramfile_ATH_TAIR10.txt 2>/dev/null';os.system(cmd)
#########################################################################################

#########################################################################################
cmd = 'time spark-submit mirLibPipeline.py ../paramfile_WHEAT_inhouse.txt 2>/dev/null';os.system(cmd)
#########################################################################################


cmd = 'free -m';os.system(cmd)




#= to study spark
#myfile = sc.textFile("file://file-path")
#myfile.saveAsTexFile("new-location")
#
#checkpoint()
#
# spark functions
#https://spark.apache.org/docs/0.6.0/api/core/spark/SparkContext.html
#https://spark.apache.org/docs/latest/configuration.html





#cmd = 'mv ../output/* ../output_depo/';os.system(cmd)


#cmd = 'free -m' #= to check memory;os.system(cmd)



#branches = 'athsettings, bowtiesplit, format_vis_input, format_vis_inputs, genomeextract, include_Dustmasker, include_RNAfold, include_bowtie, smooth, smooth_pipevarna, test, varna'.split(', ')
#for b in branches:
#  cmd = 'git branch -D ' + b
#  os.system(cmd)


#branches = 'changefoldername, faster, format_vis_input2, join, mirnauid, modifymirdup'.split(', ')
#for b in branches:
#  cmd = 'git branch -D ' + b
#  os.system(cmd)

#= download zip
#cmd = 'curl -O http://remote-server-IP/file.zip';os.system(cmd)







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


#from os import listdir
#import os.path
#rep_input = '../dbs/ATH/bowtie_split_index/'
#indirs = [f for f in listdir(rep_input) if os.path.isdir(os.path.join(rep_input, f))]
#for indir in indirs:
#  indir = rep_input + indir
#  infiles = [f for f in listdir(indir) if os.path.isfile(os.path.join(indir, f))]
#  for infile in infiles:
#    outfile = indir.rstrip('__') + '/' + infile
#    infile = indir + '/' + infile    
#    cmd = 'git mv ' + infile + ' ' + outfile
#    os.system (cmd)


#from os import listdir
#import os.path
#rep_input = '../dbs/ATH/Genome/'
#infiles = [f for f in listdir(rep_input) if os.path.isfile(os.path.join(rep_input, f))]
#for infile in infiles:
#  infile = rep_input + infile
#  outfile = infile.replace ('chr', 'Chr')
#  cmd = 'mv ' + infile + ' ' + outfile
#  os.system(cmd)



