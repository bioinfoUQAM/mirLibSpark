'''
init_dbs_customGenome.py

This script receives one file containing one to several fasta chromosomes for a genome, 
and build the bowtie index and genome folders in the format for the phase RNA pipeline.

Chao-Jung Wu
2018-11-15
'''
from __future__ import print_function
import sys, os

def makedirs_reps (reps):
  ''' mkdir in a recursive manner, meaning for a child rep, it will also make the parent reps if not exist '''
  for rep in reps:
    if not os.path.exists(rep):
      os.makedirs(rep)

def convert_multiFastaFile_to_singleFastaFiles (genomeFastaFilePath, Genome_rep):
  '''
  it reads a fasta file that contains mutiple fasta sequences, and write each of the fasta sequence a fasta file
  returns the filenames of the output single-fasta files
  '''
  with open (genomeFastaFilePath) as fh: data = [x.rstrip('\n') for x in fh.readlines()]
  flag = 0
  for i in data:
    if i.startswith('>'):
      if flag == 1: fh_out.close()
      flag = 1
      nameBase = i.split()[0][1:]
      filename = Genome_rep + nameBase + '.fa'
      fh_out = open (filename, 'w')
      print(">" + nameBase, file=fh_out)
    else: print(i, file=fh_out)

def bowtieBuild (fastaFile, NAMEbase):
  '''
  format: bowtie-build -f [path to genome fasta file] [prefix of bowtie index files]
  example: bowtie-build -f TAIR10_chr1.fas a_thaliana_t10_Chr1
  result: a_thaliana_t10_Chr1.1.ebwt and other files, total 6 files.
  '''
  cmd = 'bowtie-build -f ' + fastaFile + ' ' + NAMEbase + ' >/dev/null'
  os.system(cmd)

def create_dbs_All_bowtie (genomeFastaFilePath, rep_dbs, speciesName):
  rep_species = rep_dbs + speciesName + '/'
  bowtie_index_rep = rep_species + 'bowtie_index/'
  bowtie_index_All_rep = bowtie_index_rep + 'All/'
  reps = [bowtie_index_All_rep]
  makedirs_reps (reps)
  #
  bowtieBuild (genomeFastaFilePath, bowtie_index_All_rep + speciesName)
  cmd = 'ls ' + bowtie_index_All_rep;print('\n' + cmd);os.system(cmd)

def create_dbs_split_Genome (genomeFastaFilePath, rep_dbs, speciesName):
  rep_species = rep_dbs + speciesName + '/'
  Genome_rep = rep_species + 'Genome/'
  reps = [Genome_rep]
  makedirs_reps (reps)
  #
  convert_multiFastaFile_to_singleFastaFiles (genomeFastaFilePath, Genome_rep)
  cmd = 'ls ' + Genome_rep;print('\n' + cmd);os.system(cmd)

def init_dbs_for_a_species (genomeFastaFilePath, rep_dbs='../dbs/'):
  speciesName = genomeFastaFilePath.split("/")[-1].split('.')[0]
  create_dbs_All_bowtie (genomeFastaFilePath, rep_dbs, speciesName)
  create_dbs_split_Genome (genomeFastaFilePath, rep_dbs, speciesName)

if __name__ == '__main__' :

  if not len(sys.argv) == 2:
    sys.stderr.write('Three arguments required\nUsage: main.py [path of genome fasta file]\n')
    sys.exit()

  genomeFastaFilePath = sys.argv[1]
  init_dbs_for_a_species (genomeFastaFilePath)




