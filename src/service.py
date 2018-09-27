'''
Chao-Jung Wu
2018-09-27
'''

def createFile_KnownNonMiRNA_from_TAIR10data (infile = 'TAIR10_GFF3_genes.gff', outfile = 'TAIR10_ncRNA_CDS.gff'):
  ''' this function is not used in the pipeline, but users may use it to obtain their own KnonNonMiRNA from TAIR '''
  os.system("grep -E 'CDS|rRNA|snoRNA|snRNA|tRNA' " + infile + ' > ' + outfile)

def bowtieBuild (fastaFile, NAMEbase):
  '''
  example: bowtie-build -f TAIR10_chr1.fas a_thaliana_t10_Chr1
  result: a_thaliana_t10_Chr1.1.ebwt and other files, total 6 files.
  '''
  cmd = 'bowtie-build -f ' + fastaFile + ' ' + NAMEbase
  os.system(cmd)
