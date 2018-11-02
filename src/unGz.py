import sys
import gzip, shutil

def unzip_gzFile (gzFile):
  f_in = gzip.open(gzFile, 'rb')
  f_out = open(gzFile[:-3], 'wb')
  shutil.copyfileobj(f_in, f_out)
  f_out.close()
  f_in.close()


#gzFile = 'Triticum_aestivum.IWGSC.dna.chromosome.Un.fa.gz'
#unzip_gzFile (gzFile)

if __name__ == '__main__' :
  if not len(sys.argv) == 2:
    sys.stderr.write('main [gzFile to be unzip]')
    sys.exit()
  gzFile = sys.argv[1]
  unzip_gzFile (gzFile)


