'''
run_miranda.py

a wrapper to run miranda without spark using the function within mirLibRule.py under src folder


author: Chao-Jung Wu
date: 2018-02-12
version: 0.00.01
'''
import mirLibRules as mru

project_path='/home/cloudera/workspace/mirLibHadoop/'
miranda_exe = project_path + '/lib/miranda'

Max_Score_cutoff=170
Max_Energy_cutoff=-15
Gap_Penalty=-15

target_file = project_path + '/dbs/' + 'TAIR10_cdna_20101214_updated.fasta'
rep_tmp = project_path + '/tmp/' 

miranda_obj = mru.prog_miRanda(Max_Score_cutoff, Max_Energy_cutoff, target_file, rep_tmp, miranda_exe, Gap_Penalty)

print('test15')
