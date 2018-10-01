from os import listdir
import os.path

rep = ''
infiles = [f for f in listdir(rep) if os.path.isfile(os.path.join(rep, f))]

for infile in infiles:
  cmd = 'mv ' + infile + ' WHEAT_inhouse_' + infile
  print(cmd)
  os.system(cmd)
