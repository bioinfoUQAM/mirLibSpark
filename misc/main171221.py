def list_high_ath_mirbase ():
  High = []
  infile = 'high_conf_mature_ath_uniq100.txt'
  with open (infile, 'r') as fh:
    for i in fh:
      seq = i.split('\t')[0]
      if seq not in High:
        High.append(seq)
  print(len(High))
  return High

def list_low_ath_mirbase ():
  High = list_high_ath_mirbase ()
  Low = []
  infile = 'mirbase_all_427_real_ath.txt'
  with open (infile, 'r') as fh:
    for i in fh:
      seq = i.split('\t')[0]
      if seq not in High and seq not in Low:
        Low.append(seq)
  print(len(Low))
  return Low

def parse_12156036txt ():
  #infile = '12156036.out'
  infile = 'combined.txt'
  dictLib_mirna = {}
  with open (infile, 'r') as fh:
    #for c in range (6):
    #  fh.readline()
    line = ''
    addline = False
    for i in fh:
      if i.startswith('starting'): continue
      if i.startswith('write'): continue
      if i.startswith('line'): continue
      if i.startswith('stopping'): continue
      if i.startswith('full'): continue
      if i.startswith('  Start'): continue
      if i.startswith('NB profile_rdd'): continue

      i = i.rstrip('\n')
      if i == '  End of the processing     ': addline = False
      if addline:
        dictLib_mirna[libname].append(i)
      if i[0:28] == '--Processing of the library:':
        libname = '#' + i[30:-4]
        addline = True
        dictLib_mirna[libname] = []
  return dictLib_mirna

def uniq_mirna_inAllLibs():
  dictLib_mirna = parse_12156036txt ()
  uniqE = []
  for k, v in dictLib_mirna.items():
    for i in v:
      if i not in uniqE: uniqE.append(i)
  #print(len(uniqE))
  return uniqE

def Novel_miRNAs ():
  High = list_high_ath_mirbase ()
  Low = list_low_ath_mirbase ()
  uniqExperiment = uniq_mirna_inAllLibs()
  Novel = []
  for i in uniqExperiment:
    if not (i in High or i in Low):
      Novel.append(i)
  print(len(Novel))
  return Novel

def make_table ():
  outfile = 'tableSummaryBinary.tmp'
  fh_out = open (outfile, 'w')

  #= H = [0:100]; L = [100:349]; N = [349:1054]
  High = list_high_ath_mirbase () # 100
  Low = list_low_ath_mirbase () # 249
  Novel = sorted(Novel_miRNAs ())
  All = High + Low + Novel
  
  line = 'LIB\tnbTotal\tnbMirbaseHigh\tnbMirbaseLow\tnbNovel\t'
  for i in range(len(All)):
    if i < len(High): v = 'H_' + All[i]
    elif i < len(High+Low): v = 'L_' + All[i]
    else: v = 'N_' + All[i]
    line += v + '\t'
  line = line.rstrip('\t')
  print(line, file=fh_out, flush=True)

    
  
  libres = parse_12156036txt ()
  for k in sorted(libres.keys()):
    h, l, n = 0, 0, 0
    expressed = libres[k]

    for i in High:
      if i in expressed: h += 1
    for i in Low:
      if i in expressed: l += 1
    n = len(expressed) - h - l

    line =  k + '\t' + str(len(expressed)) + '\t' + str(h) + '\t' + str(l) + '\t' + str(n) + '\t'

    for i in All:
      if i in expressed: line += '1\t' # present
      else: line += '0\t' # absent
    line = line.rstrip('\t')
    print(line, file=fh_out, flush=True)
  fh_out.close()


      


  import transpose as tr
  infile = outfile
  outfile = infile[:-4] + '.txt'
  tr.transpose_txt(infile, outfile)
  import os
  os.remove(infile) 








#list_high_ath_mirbase ()
#list_low_ath_mirbase ()

#libres = parse_12156036txt ()
#print(len(libres))
#for k in sorted(libres.keys()):
#  print(k, len(libres[k]))

#uniq_mirna_inAllLibs()

#Novel = Novel_miRNAs ()
make_table ()
