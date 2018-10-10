'''
Chao-Jung Wu
2018-10-01

python main.py ath 1

species: 
ath: 	Arabidopsis_thaliana.TAIR10				python main.py ath 1|2
wheat: 	Triticum_aestivum.IWGSC					python main.py wheat 2
corn: 	Zea_mays.AGPv4						python main.py corn 2
rice: 	Oryza_sativa.IRGSP-1.0					python main.py rice 1
potato: Solanum_tuberosum.SolTub_3.0				python main.py potato 1
brome: 	Brachypodium_distachyon.Brachypodium_distachyon_v3.0	python main.py brome 1
wheatD: Aegilops_tauschii.ASM34733v1				python main.py wheatD 1
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

def move_chromosomeFile_to (wanted_file, ID, rep):
    fastaFile = wanted_file[:-3]
    cmd = 'rm -rf ' + wanted_file
    os.system(cmd)
    cmd = 'mv '+ fastaFile + ' ' + rep + ID + '.fa'
    os.system(cmd)

if __name__ == '__main__' :
  '''
  usage: main [organism code]
  example: python main.py ath
  '''
  if not len(sys.argv) == 3:
    sys.stderr.write('please provide organism code. \nusage: main [organism code] [1|2]\nexample: python main.py ath 1\n option 1: bowtie indexing whole genome; option 2: bowtie indexing split genome\nSupported organism codes: ath, wheat, rice, potato, brome, corn, wheatD')
    sys.exit()
  organism = sys.argv[1]
  option = sys.argv[2] #=1: All; 2: split chromosome
  
  
  if organism == 'ath':
    key = 'ATH_TAIR10'

    rep1 = '../dbs/' + key + '/'
    if not os.path.exists(rep1): os.makedirs(rep1)
    rep2 = rep1 + 'Genome/'
    if not os.path.exists(rep2): os.makedirs(rep2)
    rep3 = rep1 + 'bowtie_index/'
    if not os.path.exists(rep3): os.makedirs(rep3)

    #= get annotation files
    URL = 'ftp://ftp.ensemblgenomes.org/pub/plants/release-40/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz'
    wanted_file = 'Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz'
    curl_and_unzip_file (URL, wanted_file)
    move_file_to (wanted_file, rep1)

    URL = 'ftp://ftp.ensemblgenomes.org/pub/plants/release-40/fasta/arabidopsis_thaliana/ncrna/Arabidopsis_thaliana.TAIR10.ncrna.fa.gz'
    wanted_file = 'Arabidopsis_thaliana.TAIR10.ncrna.fa.gz'
    curl_and_unzip_file (URL, wanted_file)
    move_file_to (wanted_file, rep1)

    if option == '1':
      repAll = rep3 + 'All/'
      if not os.path.exists(repAll): os.makedirs(repAll)

      wanted_file = 'Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz'
      URL = 'ftp://ftp.ensemblgenomes.org/pub/plants/release-40/fasta/arabidopsis_thaliana/dna/' + wanted_file
      curl_and_unzip_file (URL, wanted_file)
    
      bowtieBuild (wanted_file[:-3], key)
      cmd = 'mv *.ebwt ' + repAll
      os.system(cmd)

      IDs = '1 2 3 4 5 Mt Pt'.split(' ')
      for ID in IDs:
        URL = 'ftp://ftp.ensemblgenomes.org/pub/plants/release-40/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.chromosome.'+ ID +'.fa.gz'
        wanted_file = URL.split('/')[-1]
        curl_and_unzip_file (URL, wanted_file)
        move_chromosomeFile_to (wanted_file, ID, rep2)

    if option == '2':
      IDs = '1 2 3 4 5 Mt Pt'.split(' ')
      for ID in IDs:
        rep_ch = rep3 + ID + '/'
        if not os.path.exists(rep_ch): os.makedirs(rep_ch)

        URL = 'ftp://ftp.ensemblgenomes.org/pub/plants/release-40/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.chromosome.'+ ID +'.fa.gz'
        wanted_file = URL.split('/')[-1]
        curl_and_unzip_file (URL, wanted_file)

        bowtieBuild (wanted_file[:-3], key)
        cmd = 'mv *.ebwt ' + rep_ch
        os.system(cmd)
        move_chromosomeFile_to (wanted_file, ID, rep2)
  #======================================================
  #======================================================
  #======================================================
  if organism == 'wheat':
    key = 'WHEAT_IWGSC'
    
    rep1 = '../dbs/' + key + '/'
    if not os.path.exists(rep1): os.makedirs(rep1)
    rep2 = rep1 + 'Genome/'
    if not os.path.exists(rep2): os.makedirs(rep2)
    rep3 = rep1 + 'bowtie_index/'
    if not os.path.exists(rep3): os.makedirs(rep3)

    #= get annotation files
    URL = 'ftp://ftp.ensemblgenomes.org/pub/plants/release-40/fasta/triticum_aestivum/cdna/Triticum_aestivum.IWGSC.cdna.all.fa.gz'
    wanted_file = 'Triticum_aestivum.IWGSC.cdna.all.fa.gz'
    curl_and_unzip_file (URL, wanted_file)
    move_file_to (wanted_file, rep1)

    if option == '1':
      print('option 1 not provided for wheat, switching to option 2 ...')
      option = '2'

    if option == '2':
      IDs = '1A 1B 1D 2A 2B 2D 3A 3B 3D 4A 4B 4D 5A 5B 5D 6A 6B 6D 7A 7B 7D'.split(' ')
      for ID in IDs:
        rep_ch = rep3 + ID + '/'
        if not os.path.exists(rep_ch): os.makedirs(rep_ch)

        URL = 'ftp://ftp.ensemblgenomes.org/pub/plants/release-40/fasta/triticum_aestivum/dna/Triticum_aestivum.IWGSC.dna.chromosome.'+ ID +'.fa.gz'
        wanted_file = URL.split('/')[-1]
        curl_and_unzip_file (URL, wanted_file)

        bowtieBuild (wanted_file[:-3], key)
        cmd = 'mv *.ebwt ' + rep_ch
        os.system(cmd)
        move_chromosomeFile_to (wanted_file, ID, rep2)
  #======================================================
  #======================================================
  #======================================================
  if organism == 'rice':
    key = 'RICE_IRGSP_1' 

    rep1 = '../dbs/' + key + '/'
    if not os.path.exists(rep1): os.makedirs(rep1)
    rep2 = rep1 + 'Genome/'
    if not os.path.exists(rep2): os.makedirs(rep2)
    rep3 = rep1 + 'bowtie_index/'
    if not os.path.exists(rep3): os.makedirs(rep3)

    #= get annotation files
    URL = 'ftp://ftp.ensemblgenomes.org/pub/plants/release-40/fasta/oryza_sativa/cdna/Oryza_sativa.IRGSP-1.0.cdna.all.fa.gz'
    wanted_file = URL.split('/')[-1]
    curl_and_unzip_file (URL, wanted_file)
    move_file_to (wanted_file, rep1)

    URL = 'ftp://ftp.ensemblgenomes.org/pub/plants/release-40/fasta/oryza_sativa/ncrna/Oryza_sativa.IRGSP-1.0.ncrna.fa.gz'
    wanted_file = URL.split('/')[-1]
    curl_and_unzip_file (URL, wanted_file)
    move_file_to (wanted_file, rep1)

    if option == '2':
      print('option 2 not provided for rice, switching to option 1 ...')
      option = '1'

    if option == '1':
      repAll = rep3 + 'All/'
      if not os.path.exists(repAll): os.makedirs(repAll)

      URL = 'ftp://ftp.ensemblgenomes.org/pub/plants/release-40/fasta/oryza_sativa/dna/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.gz'
      wanted_file = URL.split('/')[-1]
      curl_and_unzip_file (URL, wanted_file)
    
      bowtieBuild (wanted_file[:-3], key)
      cmd = 'mv *.ebwt ' + repAll
      os.system(cmd)

      IDs = '1 2 3 4 5 6 7 8 9 10 11 12 Mt Pt'.split(' ')
      for ID in IDs:
        URL = 'ftp://ftp.ensemblgenomes.org/pub/plants/release-40/fasta/oryza_sativa/dna/Oryza_sativa.IRGSP-1.0.dna.chromosome.'+ ID +'.fa.gz'
        wanted_file = URL.split('/')[-1]
        curl_and_unzip_file (URL, wanted_file)
        move_chromosomeFile_to (wanted_file, ID, rep2)
  #======================================================
  #======================================================
  #======================================================
  if organism == 'potato':
    key = 'POTATO_SolTub_3' 

    rep1 = '../dbs/' + key + '/'
    if not os.path.exists(rep1): os.makedirs(rep1)
    rep2 = rep1 + 'Genome/'
    if not os.path.exists(rep2): os.makedirs(rep2)
    rep3 = rep1 + 'bowtie_index/'
    if not os.path.exists(rep3): os.makedirs(rep3)

    #= get annotation files
    URL = 'ftp://ftp.ensemblgenomes.org/pub/plants/release-40/fasta/solanum_tuberosum/cdna/Solanum_tuberosum.SolTub_3.0.cdna.all.fa.gz'
    wanted_file = URL.split('/')[-1]
    curl_and_unzip_file (URL, wanted_file)
    move_file_to (wanted_file, rep1)

    URL = 'ftp://ftp.ensemblgenomes.org/pub/plants/release-40/fasta/solanum_tuberosum/ncrna/Solanum_tuberosum.SolTub_3.0.ncrna.fa.gz'
    wanted_file = URL.split('/')[-1]
    curl_and_unzip_file (URL, wanted_file)
    move_file_to (wanted_file, rep1)

    if option == '2':
      print('option 2 not provided for potato, switching to option 1 ...')
      option = '1'

    if option == '1':
      repAll = rep3 + 'All/'
      if not os.path.exists(repAll): os.makedirs(repAll)

      URL = 'ftp://ftp.ensemblgenomes.org/pub/plants/release-40/fasta/solanum_tuberosum/dna/Solanum_tuberosum.SolTub_3.0.dna.toplevel.fa.gz'
      wanted_file = URL.split('/')[-1]
      curl_and_unzip_file (URL, wanted_file)
    
      bowtieBuild (wanted_file[:-3], key)
      cmd = 'mv *.ebwt ' + repAll
      os.system(cmd)

      IDs = '1 2 3 4 5 6 7 8 9 10 11 12'.split(' ')
      for ID in IDs:
        URL = 'ftp://ftp.ensemblgenomes.org/pub/plants/release-40/fasta/solanum_tuberosum/dna/Solanum_tuberosum.SolTub_3.0.dna.chromosome.'+ ID +'.fa.gz'
        wanted_file = URL.split('/')[-1]
        curl_and_unzip_file (URL, wanted_file)
        move_chromosomeFile_to (wanted_file, ID, rep2)
      #= irregular naming of the chromosome
      URL = 'ftp://ftp.ensemblgenomes.org/pub/plants/release-40/fasta/solanum_tuberosum/dna/Solanum_tuberosum.SolTub_3.0.dna.nonchromosomal.fa.gz'
      wanted_file = URL.split('/')[-1]
      curl_and_unzip_file (URL, wanted_file)
      move_chromosomeFile_to (wanted_file, ID, rep2)
  #======================================================
  #======================================================
  #======================================================
  #= Stiff Brome
  if organism == 'brome':
    key = 'STIFFBROME' 

    rep1 = '../dbs/' + key + '/'
    if not os.path.exists(rep1): os.makedirs(rep1)
    rep2 = rep1 + 'Genome/'
    if not os.path.exists(rep2): os.makedirs(rep2)
    rep3 = rep1 + 'bowtie_index/'
    if not os.path.exists(rep3): os.makedirs(rep3)

    #= get annotation files
    URL = 'ftp://ftp.ensemblgenomes.org/pub/plants/release-40/fasta/brachypodium_distachyon/cdna/Brachypodium_distachyon.Brachypodium_distachyon_v3.0.cdna.all.fa.gz'
    wanted_file = URL.split('/')[-1]
    curl_and_unzip_file (URL, wanted_file)
    move_file_to (wanted_file, rep1)

    if option == '2':
      print('option 2 not provided for brome, switching to option 1 ...')
      option = '1'

    if option == '1':
      repAll = rep3 + 'All/'
      if not os.path.exists(repAll): os.makedirs(repAll)

      URL = 'ftp://ftp.ensemblgenomes.org/pub/plants/release-40/fasta/brachypodium_distachyon/dna/Brachypodium_distachyon.Brachypodium_distachyon_v3.0.dna.toplevel.fa.gz'
      wanted_file = URL.split('/')[-1]
      curl_and_unzip_file (URL, wanted_file)
    
      bowtieBuild (wanted_file[:-3], key)
      cmd = 'mv *.ebwt ' + repAll
      os.system(cmd)

      IDs = '1 2 3 4 5'.split(' ')
      for ID in IDs:
        URL = 'ftp://ftp.ensemblgenomes.org/pub/plants/release-40/fasta/solanum_tuberosum/dna/Solanum_tuberosum.SolTub_3.0.dna.chromosome.'+ ID +'.fa.gz'
        wanted_file = URL.split('/')[-1]
        curl_and_unzip_file (URL, wanted_file)
        move_chromosomeFile_to (wanted_file, ID, rep2)
      #= irregular naming of the chromosome
      URL = 'ftp://ftp.ensemblgenomes.org/pub/plants/release-40/fasta/brachypodium_distachyon/dna/Brachypodium_distachyon.Brachypodium_distachyon_v3.0.dna.nonchromosomal.fa.gz'
      wanted_file = URL.split('/')[-1]
      curl_and_unzip_file (URL, wanted_file)
      move_chromosomeFile_to (wanted_file, ID, rep2)
  #======================================================
  #======================================================
  #======================================================
  if organism == 'corn':
    key = 'CORN_AGPv4'
    
    rep1 = '../dbs/' + key + '/'
    if not os.path.exists(rep1): os.makedirs(rep1)
    rep2 = rep1 + 'Genome/'
    if not os.path.exists(rep2): os.makedirs(rep2)
    rep3 = rep1 + 'bowtie_index/'
    if not os.path.exists(rep3): os.makedirs(rep3)

    #= get annotation files
    URL = 'ftp://ftp.ensemblgenomes.org/pub/plants/release-40/fasta/zea_mays/cdna/Zea_mays.AGPv4.cdna.all.fa.gz'
    wanted_file = URL.split('/')[-1]
    curl_and_unzip_file (URL, wanted_file)
    move_file_to (wanted_file, rep1)

    URL = 'ftp://ftp.ensemblgenomes.org/pub/plants/release-40/fasta/zea_mays/ncrna/Zea_mays.AGPv4.ncrna.fa.gz'
    wanted_file = URL.split('/')[-1]
    curl_and_unzip_file (URL, wanted_file)
    move_file_to (wanted_file, rep1)

    if option == '1':
      print('option 1 not provided for corn, switching to option 2 ...')
      option = '2'

    if option == '2':
      IDs = '1 2 3 4 5 6 7 8 9 10 Mt Pt'.split(' ')
      for ID in IDs:
        rep_ch = rep3 + ID + '/'
        if not os.path.exists(rep_ch): os.makedirs(rep_ch)
        URL = 'ftp://ftp.ensemblgenomes.org/pub/plants/release-40/fasta/zea_mays/dna/Zea_mays.AGPv4.dna.chromosome.'+ ID +'.fa.gz'
        wanted_file = URL.split('/')[-1]
        curl_and_unzip_file (URL, wanted_file)

        bowtieBuild (wanted_file[:-3], key)
        cmd = 'mv *.ebwt ' + rep_ch
        os.system(cmd)
        move_chromosomeFile_to (wanted_file, ID, rep2)
      #= irregular naming of the chromosome
      URL = 'ftp://ftp.ensemblgenomes.org/pub/plants/release-40/fasta/zea_mays/dna/Zea_mays.AGPv4.dna.nonchromosomal.fa.gz'
      wanted_file = URL.split('/')[-1]
      curl_and_unzip_file (URL, wanted_file)
      move_chromosomeFile_to (wanted_file, ID, rep2)
  #======================================================
  #======================================================
  #======================================================
  if organism == 'wheatD':
    key = 'WHEAT_D'
    
    rep1 = '../dbs/' + key + '/'
    if not os.path.exists(rep1): os.makedirs(rep1)
    rep2 = rep1 + 'Genome/'
    if not os.path.exists(rep2): os.makedirs(rep2)
    rep3 = rep1 + 'bowtie_index/'
    if not os.path.exists(rep3): os.makedirs(rep3)

    #= get annotation files
    URL = 'ftp://ftp.ensemblgenomes.org/pub/plants/release-40/fasta/aegilops_tauschii/cdna/Aegilops_tauschii.ASM34733v1.cdna.all.fa.gz'
    wanted_file = URL.split('/')[-1]
    curl_and_unzip_file (URL, wanted_file)
    move_file_to (wanted_file, rep1)

    URL = 'ftp://ftp.ensemblgenomes.org/pub/plants/release-40/fasta/aegilops_tauschii/ncrna/Aegilops_tauschii.ASM34733v1.ncrna.fa.gz'
    wanted_file = URL.split('/')[-1]
    curl_and_unzip_file (URL, wanted_file)
    move_file_to (wanted_file, rep1)

    if option == '2':
      print('option 2 not provided for wheatD, switching to option 1 ...')
      option = '1'

    if option == '1':
      repAll = rep3 + 'All/'
      if not os.path.exists(repAll): os.makedirs(repAll)

      URL = 'ftp://ftp.ensemblgenomes.org/pub/plants/release-40/fasta/aegilops_tauschii/dna/Aegilops_tauschii.ASM34733v1.dna.toplevel.fa.gz'
      wanted_file = URL.split('/')[-1]
      curl_and_unzip_file (URL, wanted_file)
    
      bowtieBuild (wanted_file[:-3], key)
      cmd = 'mv *.ebwt ' + repAll
      os.system(cmd)

      #= irregular naming of the chromosome
      URL = 'ftp://ftp.ensemblgenomes.org/pub/plants/release-40/fasta/aegilops_tauschii/dna/Aegilops_tauschii.ASM34733v1.dna.toplevel.fa.gz'
      wanted_file = URL.split('/')[-1]
      curl_and_unzip_file (URL, wanted_file)
      move_chromosomeFile_to (wanted_file, ID, rep2)
      #
      URL = 'ftp://ftp.ensemblgenomes.org/pub/plants/release-40/fasta/aegilops_tauschii/dna/Aegilops_tauschii.ASM34733v1.dna.nonchromosomal.fa.gz'
      wanted_file = URL.split('/')[-1]
      curl_and_unzip_file (URL, wanted_file)
      move_chromosomeFile_to (wanted_file, ID, rep2)

