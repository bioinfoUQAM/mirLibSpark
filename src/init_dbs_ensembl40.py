'''
Chao-Jung Wu
2018-10-01

option1: small genomes (< 1G), or the genomes that ensembl does not provide chromosomal files
option2: large genomes (> 1G)
== supported species ========================================== cmd line ============================== test ================
ath: 	Arabidopsis_thaliana.TAIR10				python init_dbs_ensembl40.py ath 1|2	tested_previously_ok
wheat: 	Triticum_aestivum.IWGSC					python init_dbs_ensembl40.py wheat 2	tested_previously_ok
corn: 	Zea_mays.AGPv4						python init_dbs_ensembl40.py corn 2
rice: 	Oryza_sativa.IRGSP-1.0					python init_dbs_ensembl40.py rice 1	tested_181010_ok
potato: Solanum_tuberosum.SolTub_3.0				python init_dbs_ensembl40.py potato 1   tested_181010_ok
brome: 	Brachypodium_distachyon.Brachypodium_distachyon_v3.0	python init_dbs_ensembl40.py brome 1    tested_181011_ok
wheatD: Aegilops_tauschii.ASM34733v1				python init_dbs_ensembl40.py wheatD 1
'''


import sys, os
import gzip, shutil




ensembl_version = 'release-40'





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
  cmd = 'bowtie-build -f ' + fastaFile + ' ' + NAMEbase
  os.system(cmd)

def curl_and_unzip_file (place, URL, wanted_file):
    if place == 'S1' or place == 'L':
      cmd = 'curl -O ' + URL
      os.system(cmd)
    unzip_gzFile (wanted_file)

def move_file_to (wanted_file, rep):
    fastaFile = wanted_file[:-3]
    cmd = 'rm -rf ' + wanted_file
    os.system(cmd)
    cmd = 'mv '+ fastaFile + ' ' + rep
    os.system(cmd)

def move_chromosomeFile_to (wanted_file, ID, rep):
    fastaFile = wanted_file[:-3]
    cmd = 'rm -rf ' + wanted_file
    os.system(cmd)
    cmd = 'mv '+ fastaFile + ' ' + rep + ID + '.fa'
    os.system(cmd)

def option2_forLargeGenome (place, organism, key, URL_cdna, URL_ncrna, IDs, URL_dna_prefix, URL_irregulars):
    ''' generate chromosome-split genome and bowtie index '''
    rep1 = '../dbs/' + key + '/'
    if not os.path.exists(rep1): os.makedirs(rep1)
    rep2 = rep1 + 'Genome/'
    if not os.path.exists(rep2): os.makedirs(rep2)
    rep3 = rep1 + 'bowtie_index/'
    if not os.path.exists(rep3): os.makedirs(rep3)

    #= get annotation files: cdna, ncrna
    wanted_file = URL_cdna.split('/')[-1]
    curl_and_unzip_file (place, URL_cdna, wanted_file)
    move_file_to (wanted_file, rep1)

    if not URL_ncrna == '':
      wanted_file = URL_ncrna.split('/')[-1]
      curl_and_unzip_file (place, URL_ncrna, wanted_file)
      move_file_to (wanted_file, rep1)

    #= build chromosome-split bowtie index in chromosome-specific sub-folders
    #= get chromosome-split genome in Genome folders
    for ID in IDs:
      rep_ch = rep3 + ID + '/'
      if not os.path.exists(rep_ch): os.makedirs(rep_ch)

      URL_dna = URL_dna_prefix + ID +'.fa.gz'
      wanted_file = URL_dna.split('/')[-1]
      curl_and_unzip_file (place, URL_dna, wanted_file)

      #''' #skip build
      if place == 'L' or place == 'S2':
        bowtieBuild (wanted_file[:-3], key)
        cmd = 'mv *.ebwt ' + rep_ch
        os.system(cmd)
        move_chromosomeFile_to (wanted_file, ID, rep2)
      #'''

    #= irregular naming of the chromosomes
    for URL in URL_irregulars:
      wanted_file = URL.split('/')[-1]
      curl_and_unzip_file (place, URL, wanted_file)
      ID = wanted_file.split('.')[0]
      move_chromosomeFile_to (wanted_file, ID, rep2)

def option1_forSmallGenome (place, organism, key, URL_cdna, URL_ncrna, URL_toplevel, IDs, URL_dna_prefix, URL_irregulars):
    ''' generate whole-genome bowtie index and split chromosome genome '''
    rep1 = '../dbs/' + key + '/'
    if not os.path.exists(rep1): os.makedirs(rep1)
    rep2 = rep1 + 'Genome/'
    if not os.path.exists(rep2): os.makedirs(rep2)
    rep3 = rep1 + 'bowtie_index/'
    if not os.path.exists(rep3): os.makedirs(rep3)

    #= get annotation files: cdna, ncrna
    wanted_file = URL_cdna.split('/')[-1]
    curl_and_unzip_file (place, URL_cdna, wanted_file)
    move_file_to (wanted_file, rep1)

    if not URL_ncrna == '':
      wanted_file = URL_ncrna.split('/')[-1]
      curl_and_unzip_file (place, URL_ncrna, wanted_file)
      move_file_to (wanted_file, rep1)

    repAll = rep3 + 'All/'
    if not os.path.exists(repAll): os.makedirs(repAll)

    wanted_file = URL_toplevel.split('/')[-1]
    curl_and_unzip_file (place, URL_toplevel, wanted_file)
    
    #''' #skip build
    if place == 'L' or place == 'S2':
      bowtieBuild (wanted_file[:-3], key)
      cmd = 'mv *.ebwt ' + repAll
      os.system(cmd)
    #'''

    for ID in IDs:
      URL_dna = URL_dna_prefix + ID +'.fa.gz'
      wanted_file = URL_dna.split('/')[-1]
      curl_and_unzip_file (place, URL_dna, wanted_file)
      move_chromosomeFile_to (wanted_file, ID, rep2)

    #= irregular naming of the chromosomes
    for URL in URL_irregulars:
      wanted_file = URL.split('/')[-1]
      curl_and_unzip_file (place, URL, wanted_file)
      ID = wanted_file.split('.')[0]
      move_chromosomeFile_to (wanted_file, ID, rep2)

#=========================
#
#  MAIN ()
#
#=========================
if __name__ == '__main__' :
  '''
  usage: init_dbs_ensembl40 [organism code] [1|2] [L|S1|S2]
  example: python init_dbs_ensembl40.py ath 1
  '''
  if not len(sys.argv) == 4:
    sys.stderr.write('usage: main [organism code] [1|2] [L|S1|S2]\nexample: python init_dbs_ensembl40.py ath 1 L\n option 1: bowtie indexing whole genome; option 2: bowtie indexing split genome\nSupported organism codes: ath, wheat, rice, potato, brome, corn, wheatD')
    sys.exit()
  organism = sys.argv[1]
  option = sys.argv[2] #=1: All; 2: split chromosome
  place = sys.argv[3]  #= L: local (do curl and bowtie-build); S1: server, do only curl, not bowtie-build; S2: server, do only bowtie-build, not curl
  
  if organism == 'ath':
    key = 'ATH_TAIR10'
    URL_cdna = 'ftp://ftp.ensemblgenomes.org/pub/plants/' + ensembl_version + '/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz'
    URL_ncrna = 'ftp://ftp.ensemblgenomes.org/pub/plants/' + ensembl_version + '/fasta/arabidopsis_thaliana/ncrna/Arabidopsis_thaliana.TAIR10.ncrna.fa.gz'
    URL_toplevel = 'ftp://ftp.ensemblgenomes.org/pub/plants/' + ensembl_version + '/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz' 
    IDs =  '1 2 3 4 5 Mt Pt'.split(' ')
    URL_dna_prefix = 'ftp://ftp.ensemblgenomes.org/pub/plants/' + ensembl_version + '/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.chromosome.'
    URL_irregulars = []
    if option == '1': option1_forSmallGenome (place, organism, key, URL_cdna, URL_ncrna, URL_toplevel, IDs, URL_dna_prefix, URL_irregulars)
    elif option == '2': option2_forLargeGenome (place, organism, key, URL_cdna, URL_ncrna, IDs, URL_dna_prefix, URL_irregulars)
  #======================================================
  if organism == 'wheat':
    key = 'WHEAT_IWGSC'
    URL_cdna = 'ftp://ftp.ensemblgenomes.org/pub/plants/' + ensembl_version + '/fasta/triticum_aestivum/cdna/Triticum_aestivum.IWGSC.cdna.all.fa.gz'
    URL_ncrna = ''
    IDs = '1A 1B 1D 2A 2B 2D 3A 3B 3D 4A 4B 4D 5A 5B 5D 6A 6B 6D 7A 7B 7D'.split(' ')
    URL_dna_prefix = 'ftp://ftp.ensemblgenomes.org/pub/plants/' + ensembl_version + '/fasta/triticum_aestivum/dna/Triticum_aestivum.IWGSC.dna.chromosome.'
    URL_irregulars = []
    if option == '1':
      print('option 1 not provided for ' + organism + ', switching to option 2 ...')
      option = '2'
    option2_forLargeGenome (place, organism, key, URL_cdna, URL_ncrna, IDs, URL_dna_prefix, URL_irregulars)
  #======================================================
  if organism == 'corn':
    key = 'CORN_AGPv4'
    URL_cdna = 'ftp://ftp.ensemblgenomes.org/pub/plants/' + ensembl_version + '/fasta/zea_mays/cdna/Zea_mays.AGPv4.cdna.all.fa.gz'
    URL_ncrna = 'ftp://ftp.ensemblgenomes.org/pub/plants/' + ensembl_version + '/fasta/zea_mays/ncrna/Zea_mays.AGPv4.ncrna.fa.gz'
    IDs = '1 2 3 4 5 6 7 8 9 10 Mt Pt'.split(' ')
    URL_dna_prefix = 'ftp://ftp.ensemblgenomes.org/pub/plants/' + ensembl_version + '/fasta/zea_mays/dna/Zea_mays.AGPv4.dna.chromosome.'
    URL_irregulars = ['ftp://ftp.ensemblgenomes.org/pub/plants/' + ensembl_version + '/fasta/zea_mays/dna/Zea_mays.AGPv4.dna.nonchromosomal.fa.gz']
    if option == '1':
      print('option 1 not provided for ' + organism + ', switching to option 2 ...')
      option = '2'
    option2_forLargeGenome (place, organism, key, URL_cdna, URL_ncrna, IDs, URL_dna_prefix, URL_irregulars)
  #======================================================
  if organism == 'rice':
    key = 'RICE_IRGSP_1' 
    URL_cdna = 'ftp://ftp.ensemblgenomes.org/pub/plants/' + ensembl_version + '/fasta/oryza_sativa/cdna/Oryza_sativa.IRGSP-1.0.cdna.all.fa.gz'
    URL_ncrna = 'ftp://ftp.ensemblgenomes.org/pub/plants/' + ensembl_version + '/fasta/oryza_sativa/ncrna/Oryza_sativa.IRGSP-1.0.ncrna.fa.gz'
    URL_toplevel = 'ftp://ftp.ensemblgenomes.org/pub/plants/' + ensembl_version + '/fasta/oryza_sativa/dna/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.gz'
    IDs = '1 2 3 4 5 6 7 8 9 10 11 12 Mt Pt'.split(' ')
    URL_dna_prefix = 'ftp://ftp.ensemblgenomes.org/pub/plants/' + ensembl_version + '/fasta/oryza_sativa/dna/Oryza_sativa.IRGSP-1.0.dna.chromosome.'
    URL_irregulars = []
    if option == '2':
      print('option 2 not provided for ' + organism + ', switching to option 1 ...')
      option = '1'
    option1_forSmallGenome (place, organism, key, URL_cdna, URL_ncrna, URL_toplevel, IDs, URL_dna_prefix, URL_irregulars)
  #======================================================
  if organism == 'potato':
    key = 'POTATO_SolTub_3' 
    URL_cdna = 'ftp://ftp.ensemblgenomes.org/pub/plants/' + ensembl_version + '/fasta/solanum_tuberosum/cdna/Solanum_tuberosum.SolTub_3.0.cdna.all.fa.gz'
    URL_ncrna = 'ftp://ftp.ensemblgenomes.org/pub/plants/' + ensembl_version + '/fasta/solanum_tuberosum/ncrna/Solanum_tuberosum.SolTub_3.0.ncrna.fa.gz'
    URL_toplevel = 'ftp://ftp.ensemblgenomes.org/pub/plants/' + ensembl_version + '/fasta/solanum_tuberosum/dna/Solanum_tuberosum.SolTub_3.0.dna.toplevel.fa.gz'
    IDs = '1 2 3 4 5 6 7 8 9 10 11 12'.split(' ')   
    URL_dna_prefix = 'ftp://ftp.ensemblgenomes.org/pub/plants/' + ensembl_version + '/fasta/solanum_tuberosum/dna/Solanum_tuberosum.SolTub_3.0.dna.chromosome.'
    URL_irregulars = ['ftp://ftp.ensemblgenomes.org/pub/plants/' + ensembl_version + '/fasta/solanum_tuberosum/dna/Solanum_tuberosum.SolTub_3.0.dna.nonchromosomal.fa.gz']
    if option == '2':
      print('option 2 not provided for ' + organism + ', switching to option 1 ...')
      option = '1'
    option1_forSmallGenome (place, organism, key, URL_cdna, URL_ncrna, URL_toplevel, IDs, URL_dna_prefix, URL_irregulars)
  #======================================================
  if organism == 'brome':
    key = 'STIFFBROME' 
    URL_cdna = 'ftp://ftp.ensemblgenomes.org/pub/plants/' + ensembl_version + '/fasta/brachypodium_distachyon/cdna/Brachypodium_distachyon.Brachypodium_distachyon_v3.0.cdna.all.fa.gz'
    URL_ncrna = ''
    URL_toplevel = 'ftp://ftp.ensemblgenomes.org/pub/plants/' + ensembl_version + '/fasta/brachypodium_distachyon/dna/Brachypodium_distachyon.Brachypodium_distachyon_v3.0.dna.toplevel.fa.gz'
    IDs = IDs = '1 2 3 4 5'.split(' ')   
    URL_dna_prefix = 'ftp://ftp.ensemblgenomes.org/pub/plants/' + ensembl_version + '/fasta/brachypodium_distachyon/dna/Brachypodium_distachyon.Brachypodium_distachyon_v3.0.dna.chromosome.'
    URL_irregulars = ['ftp://ftp.ensemblgenomes.org/pub/plants/' + ensembl_version + '/fasta/brachypodium_distachyon/dna/Brachypodium_distachyon.Brachypodium_distachyon_v3.0.dna.nonchromosomal.fa.gz']
    if option == '2':
      print('option 2 not provided for ' + organism + ', switching to option 1 ...')
      option = '1'
    option1_forSmallGenome (place, organism, key, URL_cdna, URL_ncrna, URL_toplevel, IDs, URL_dna_prefix, URL_irregulars)
  #======================================================
  if organism == 'wheatD':
    key = 'WHEAT_D'
    URL_cdna = 'ftp://ftp.ensemblgenomes.org/pub/plants/' + ensembl_version + '/fasta/aegilops_tauschii/cdna/Aegilops_tauschii.ASM34733v1.cdna.all.fa.gz'    
    URL_ncrna = 'ftp://ftp.ensemblgenomes.org/pub/plants/' + ensembl_version + '/fasta/aegilops_tauschii/ncrna/Aegilops_tauschii.ASM34733v1.ncrna.fa.gz'
    URL_toplevel = 'ftp://ftp.ensemblgenomes.org/pub/plants/' + ensembl_version + '/fasta/aegilops_tauschii/dna/Aegilops_tauschii.ASM34733v1.dna.toplevel.fa.gz'
    IDs = []  
    URL_dna_prefix = ''
    URL_irregulars = ['ftp://ftp.ensemblgenomes.org/pub/plants/' + ensembl_version + '/fasta/aegilops_tauschii/dna/Aegilops_tauschii.ASM34733v1.dna.toplevel.fa.gz', \
                      'ftp://ftp.ensemblgenomes.org/pub/plants/' + ensembl_version + '/fasta/aegilops_tauschii/dna/Aegilops_tauschii.ASM34733v1.dna.nonchromosomal.fa.gz']
    if option == '2':
      print('option 2 not provided for ' + organism + ', switching to option 1 ...')
      option = '1'
    option1_forSmallGenome (place, organism, key, URL_cdna, URL_ncrna, URL_toplevel, IDs, URL_dna_prefix, URL_irregulars)

