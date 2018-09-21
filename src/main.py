import os

#cmd = 'python purge.py'
#os.system (cmd)

#cmd = 'git pull origin bowtiesplit'
#os.system (cmd)

#cmd = 'rm -fr ../input/fake_a3.txt'
#os.system(cmd)


#cmd = 'cp ../input_samples/fake_a3.txt ../input'
#os.system(cmd)




cmd = 'time spark-submit mirLibPipeline.py ../paramfile.txt 2>/dev/null'
os.system(cmd)


#cmd = 'mkdir ../dbs/ATH/bowtie_split_index'
#os.system(cmd)


#chro = [1, 2, 3, 4, 5, 'C', 'M']
#for c in chro:
#  refgenome = '../dbs/ATH/Genome/TAIR10_chr' + str(c) + '.fas'
#  cmd = '../lib/bowtie-build -f ' + refgenome + ' atht10_' + 'chr' + str(c)
#  os.system(cmd)


