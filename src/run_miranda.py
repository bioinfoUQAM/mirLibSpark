'''
run_miranda.py

a standalone wrapper to run miranda without spark


author: Chao-Jung Wu
date: 2018-02-12
version: 0.00.01
'''
import os
from operator import itemgetter





class prog_miRanda ():
  def __init__ (self, Max_Score_cutoff, Max_Energy_cutoff, target_file, rep_tmp, miranda_exe, Gap_Penalty):
    #self.env = os.environ

    #== variables ==
    self.Max_Score_cutoff = Max_Score_cutoff
    self.Max_Energy_cutoff = Max_Energy_cutoff
    self.Gap_Penalty = Gap_Penalty

    self.target_file = target_file
    self.rep_tmp = rep_tmp
    self.miranda_exe = miranda_exe

    self.dict_seq_target = {}

  def computeTargetbyMiranda (self, e):
    '''
    $miranda examples/bantam_stRNA.fasta examples/hid_UTR.fasta
    $miranda ../tmp/tmp_mirna_seq.txt ../Arabidopsis/TAIR/Genome/TAIR10_blastsets/TAIR10_cdna_20101214_updated_1cdna.fasta
    $miranda ../tmp/GCTCACTGCTCTTTCTGTCAGA_tmpseq_forMiranda.txt TAIR10_cdna_20101214_updated_39cdna.fasta

    ## NOTE before disable miranda (170714): need to modify the code to use options such as -sc, -en, -go, -ge, -quiet
    '''

    miRNAseq = e #= be careful "return e" and "e[1].append(target_results)"

    '''
    if e[0] in self.dict_seq_target.keys():
      e[1].append(self.dict_seq_target[e[0]])
      return e
    #'''

    tmp_file = self.rep_tmp + miRNAseq + '_tmpseq_forMiranda2.txt' 
    with open (tmp_file, 'w') as fh_tmp: print >> fh_tmp, '>x\n' + miRNAseq

    #cmd = [self.miranda_exe, tmp_file, self.target_file, '-strict', '-sc', self.Max_Score_cutoff, '-en', self.Max_Energy_cutoff, '-go', self.Gap_Penalty]

    cmd = self.miranda_exe + ' ' + tmp_file + ' ' + self.target_file + ' -strict -sc ' + str(self.Max_Score_cutoff) + ' -en ' + str(self.Max_Energy_cutoff) + ' -go ' + str(self.Gap_Penalty) + '> tmpmiranda2.txt'
    os.system(cmd)

    with open ('tmpmiranda.txt', 'r') as fh:
      mirandaout = fh.readlines()

    target_results = []

    for i in mirandaout[30:]: 
    #= because the first 30ish lines contain only program description
      if i[:3] == '>>x': 
        #= isTargteet == [Seq1, Seq2, Tot_Score, Tot_Energy, Max_Score, Max_Energy, Strand, Len1, Len2, Positions]
        target_result = i.split('\t') 
        target_results.append(target_result[1:])
        #target_results.append([target_result[1], target_result[9]]) #= only record gene and positions

    #= target_results == [[target1], [target2], ...]
    #= [['AT1G51370.2', '306.00', '-36.41', '153.00', '-20.70', '1', '23', '1118', ' 20 698']]
    #= [gene, total_score, total_energy, max_score, max_energy, strand, len_miRNA, len_gene, postions]

    #= targets are sorted from highest total scores to lowest total scores
    target_results = sorted(target_results, key=itemgetter(1), reverse=True)
    #= only the top 15 targets are curated for report
    if len(target_results) > 15: target_results = target_results[:15]
    self.dict_seq_target[miRNAseq] = target_results
    #e[1].append(target_results)
    #return e
    return target_results

def parse():
  infile = '../input/predicted859mirna.txt'
  master_predicted_mirna = []
  with open (infile, 'r') as fh:
    for i in fh:
      mirna = i.rstrip('\n')
      if mirna not in master_predicted_mirna: master_predicted_mirna.append(mirna)
  return master_predicted_mirna


project_path='/home/cloudera/workspace/mirLibHadoop/'
miranda_exe = project_path + '/lib/miranda'

Max_Score_cutoff=170
Max_Energy_cutoff=-15
Gap_Penalty=-15

target_file = project_path + '/dbs/' + 'TAIR10_cdna_20101214_updated.fasta'
rep_tmp = project_path + '/tmp/' 

miranda_obj = prog_miRanda(Max_Score_cutoff, Max_Energy_cutoff, target_file, rep_tmp, miranda_exe, Gap_Penalty)

outfile = 'tgtg.temp.txt'
fh_out = open (outfile, 'w')

#miRNA = 'TCATGGTCAGATCCGTCATCC'
master_predicted_mirna = parse()
for miRNA in master_predicted_mirna:
  print >> fh_out, '>' + miRNA
  target_results = miranda_obj.computeTargetbyMiranda(miRNA)
  for i in target_results: print >> fh_out, i
  break

fh_out.close()

print('test15')
