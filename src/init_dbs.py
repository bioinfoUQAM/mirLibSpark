'''
Chao-Jung Wu
2018-10-01

python main.py ath

'''


import sys, os
import gzip, shutil

def unzip_gzFile (gzFile):
  f_in = gzip.open(gzFile, 'rb')
  f_out = open(gzFile[:-3], 'wb')
  shutil.copyfileobj(f_in, f_out)
  f_out.close()
  f_in.close()

def bowtieBuild (fastaFile, NAMEbase):
  '''
  example: bowtie-build -f TAIR10_chr1.fas a_thaliana_t10_Chr1
  result: a_thaliana_t10_Chr1.1.ebwt and other files, total 6 files.
  '''
  cmd = '../lib/bowtie-build -f ' + fastaFile + ' ' + NAMEbase
  os.system(cmd)

def curl_and_unzip_file (URL, wanted_file):
    cmd = 'curl -O ' + URL
    os.system(cmd)
    unzip_gzFile (wanted_file)

def move_file_to (wanted_file, rep):
    fastaFile = wanted_file[:-3]
    cmd = 'rm -rf ' + wanted_file
    os.system(cmd)
    cmd = 'mv '+ fastaFile + ' ' + rep
    os.system(cmd)

if __name__ == '__main__' :
  '''
  usage: main [organism 3-letter code]
  example: python main.py ath
  '''
  if not len(sys.argv) == 3:
    sys.stderr.write('please provide organism 3-letter code. \nusage: main [organism 3-letter code] [option]\nexample: python main.py ath 1')
    sys.exit()
  organism = sys.argv[1]
  option = sys.argv[2] #=1: All; 2: split chromosome
  
  
  if organism == 'ath' and option == '1':
    key = 'ATH_TAIR10'

    wanted_file = 'Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz'
    URL = 'ftp://ftp.ensemblgenomes.org/pub/plants/release-40/fasta/arabidopsis_thaliana/dna/' + wanted_file

    rep1 = '../dbs/' + key + '/'
    if not os.path.exists(rep1): os.makedirs(rep1)
    rep2 = rep1 + 'Genome/'
    if not os.path.exists(rep2): os.makedirs(rep2)
    rep3 = rep1 + 'bowtie_index/'
    if not os.path.exists(rep3): os.makedirs(rep3)
    rep3 = rep3 + 'All/'
    if not os.path.exists(rep3): os.makedirs(rep3)

    curl_and_unzip_file (URL, wanted_file)
    
    #bowtieBuild (wanted_file[:-3], key)
    cmd = 'mv *.ebwt ' + rep3
    os.system(cmd)
    #move_file_to (wanted_file, rep2)

    #= get annotation files
    URL = 'ftp://ftp.ensemblgenomes.org/pub/plants/release-40/fasta/arabidopsis_thaliana/cds/Arabidopsis_thaliana.TAIR10.cds.all.fa.gz'
    wanted_file = 'Arabidopsis_thaliana.TAIR10.cds.all.fa.gz'
    #curl_and_unzip_file (URL, wanted_file)
    #move_file_to (wanted_file, rep1)

    URL = 'ftp://ftp.ensemblgenomes.org/pub/plants/release-40/fasta/arabidopsis_thaliana/ncrna/Arabidopsis_thaliana.TAIR10.ncrna.fa.gz'
    wanted_file = 'Arabidopsis_thaliana.TAIR10.ncrna.fa.gz'
    #curl_and_unzip_file (URL, wanted_file)
    #move_file_to (wanted_file, rep1)

  ###############################################################################
  if organism == 'ath' and option == '2':
    key = 'ATH_TAIR10'
    IDs = '1 2 3 4 5 Mt Pt'.split(' ')

    rep1 = '../dbs/' + key + '/'
    if not os.path.exists(rep1): os.makedirs(rep1)
    rep2 = rep1 + 'Genome/'
    if not os.path.exists(rep2): os.makedirs(rep2)
    rep3 = rep1 + 'bowtie_index/'
    if not os.path.exists(rep3): os.makedirs(rep3)

    for ID in IDs:
      URL = 'ftp://ftp.ensemblgenomes.org/pub/plants/release-40/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.chromosome.'+ ID +'.fa.gz'
      wanted_file = URL.split('/')[-1]
      curl_and_unzip_file (URL, wanted_file)

      rep_ch = rep3 + ID + '/'
      if not os.path.exists(rep_ch): os.makedirs(rep_ch)

      bowtieBuild (wanted_file[:-3], key)
      cmd = 'mv *.ebwt ' + rep_ch
      os.system(cmd)
      move_file_to (wanted_file, rep2)

    #= get annotation files
    URL = 'ftp://ftp.ensemblgenomes.org/pub/plants/release-40/fasta/arabidopsis_thaliana/cds/Arabidopsis_thaliana.TAIR10.cds.all.fa.gz'
    wanted_file = 'Arabidopsis_thaliana.TAIR10.cds.all.fa.gz'
    curl_and_unzip_file (URL, wanted_file)
    move_file_to (wanted_file, rep1)

    URL = 'ftp://ftp.ensemblgenomes.org/pub/plants/release-40/fasta/arabidopsis_thaliana/ncrna/Arabidopsis_thaliana.TAIR10.ncrna.fa.gz'
    wanted_file = 'Arabidopsis_thaliana.TAIR10.ncrna.fa.gz'
    curl_and_unzip_file (URL, wanted_file)
    move_file_to (wanted_file, rep1)


