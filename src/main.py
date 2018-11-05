'''
Graham manual 2018-10-26

modify `$paramfile_ATH_TAIR10_graham.txt`
do `$mkdir input` in project folder,
do `$cp input_samples/100.txt input/`
do `$ cd workdir`,
do `$mkdir jobout` in workdir,
do `$sbatch pyspark_submit_jv2_20_ath032.sh`
'''

import os

#= use branch "jmajor" as my own master to develop features for managing the memory

#cmd = 'sudo mount -t vboxsf Projects /home/cloudera/vm_dropbox_projects/';os.system(cmd)
#cmd = 'sudo mount -t vboxsf VM_Share /home/cloudera/vm_share/';os.system(cmd)


#cmd = 'python purge.py';os.system (cmd)
#cmd = 'git pull origin graham';os.system (cmd)
#cmd = 'git branch -D dbsauto';os.system(cmd)
#cmd = 'git remote prune origin';os.system(cmd)
#cmd = 'git fetch origin';os.system(cmd)
#cmd = 'git checkout master';os.system(cmd)


#=== ATH diff =======
cmd = 'cp ../input_samples/fake_a.txt ../input';os.system(cmd)
#cmd = 'rm -fr ../input/fake_a.txt';os.system(cmd)

cmd = 'cp ../input_samples/fake_a3.txt ../input';os.system(cmd)
#cmd = 'rm -fr ../input/fake_a3.txt';os.system(cmd)
#====================


#=== WHEAT diff =======
#cmd = 'cp ../input_samples/wheat_3A_a1.txt ../input';os.system(cmd)
#cmd = 'rm -fr ../input/wheat_3A_a1.txt';os.system(cmd)

#cmd = 'cp ../input_samples/wheat_3A_a2.txt ../input';os.system(cmd)
#cmd = 'rm -fr ../input/wheat_3A_a2.txt';os.system(cmd)
#====================



#cmd = 'cp ../input_samples/high_conf_mature_ath_uniq_raw.txt ../input';os.system(cmd)
#cmd = 'rm -fr ../input/high_conf_mature_ath_uniq_raw.txt';os.system(cmd)


#cmd = 'cp ../input_samples/fake_a5.txt ../input';os.system(cmd)
#cmd = 'rm -fr ../input/fake_a5.txt';os.system(cmd)

#cmd = 'cp ../input_samples/100.txt ../input';os.system(cmd)
#cmd = 'rm -f ../input/100.txt';os.system(cmd)

#cmd = 'cp ../input_samples/11w2013_t2_1.fasta ../input';os.system(cmd)
#cmd = 'rm ../input/11w2013_t2_1.fasta ../input_samples';os.system(cmd)



#cmd = 'git add -f ../output/*';os.system(cmd)
#cmd = 'git commit -m "save a copy of output"';os.system(cmd)
#cmd = 'git push origin refactor';os.system(cmd)



#cmd = 'free -m';os.system(cmd)

#########################################################################################
cmd = 'time spark-submit mirLibPipeline.py 2>/dev/null';os.system(cmd)
#########################################################################################

#########################################################################################
#cmd = 'time spark-submit mirLibPipeline.py --species wheat 2>/dev/null';os.system(cmd)
#########################################################################################


#cmd = 'free -m';os.system(cmd)




#sed -n 2000,2004p file

#dos2unix windows.txt unix.txt


#= to study spark
#myfile = sc.textFile("file://file-path")
#myfile.saveAsTexFile("new-location")
#
#checkpoint()
#
# spark functions
#https://spark.apache.org/docs/0.6.0/api/core/spark/SparkContext.html
#https://spark.apache.org/docs/latest/configuration.html


##Remove CTRL-M characters from a file in UNIX Inside vi [in ESC mode] type: 
#:%s/^M//g
#To enter ^M, type CTRL-V, then CTRL-M. That is, hold down the CTRL key then press V and M in succession.


#cmd = 'mv ../output/* ../output_depo/';os.system(cmd)


#cmd = 'bowtie-build -f XLOC_126981_6D_REVCOMP.fasta XLOC_126981_6D_REVCOMP';os.system(cmd)


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

#path note: 181002
#/home/cloudera/vm_dropbox_projects/1_Project_mirLibHadoop/gitRepo_mirLibHadoop/
#/home/cloudera/workspace/


