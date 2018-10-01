'''
Chao-Jung Wu
2018-09-26

#= ats: Aegilops tauschii (wheat D)
#= sot: potato, Solanum tuberosum

'''

import sys, os

#organism = 'ath, osa, zma'.split(', ')

if __name__ == '__main__' :
  '''
  usage: main [organim 3-letter code]
  example: python main.py ath
  '''
  if not len(sys.argv) == 2:
    sys.stderr.write('please provide organim 3-letter code. \nusage: main [organim 3-letter code]\nexample: python main.py ath')
    sys.exit()
  organism = sys.argv[1]

  if not os.path.exists('../output/'): os.makedirs('../output/')
  cmd = 'python bioservicesKEGG.py ' + organism
  os.system(cmd)

  rep = '../output/' + organism + '/'
  outfile = '../output/' + organism + '_gene_vs_pathway.txt'
  outfile2 = '../output/' + organism + '_pathway_description.txt'
  outfile3 = '../output/' + organism + '_gene_description.txt'
  fh_out = open (outfile, 'w')
  fh_out2 = open (outfile2, 'w')
  fh_out3 = open (outfile3, 'w')
  infiles  = [rep + f for f in os.listdir(rep) if os.path.isfile(os.path.join(rep, f))]
  for infile in infiles:
    if infile.endswith('.swp'): continue
    pathwayname = infile.split('/')[-1][:-4] 
    turnon = 0
    count = 0
    with open (infile, 'r') as fh:
      for line in fh:
        count += 1
        if count == 2: 
          desc = line[12:].split(' - ')[0]
          aline = '\t'.join( [pathwayname, desc] )
          print(aline, file=fh_out2, flush=True)
        if line.startswith('GENE'): turnon = 1
        if line.startswith('COMPOUND'): turnon = 0
        if turnon == 1:
          data = line[12:].rstrip('\n').split('  ')
          if not len(data) == 2:continue 
          geneid = data[0]
          geneinfo = data[1]
          line = '\t'.join([geneid, pathwayname])
          print(line, file=fh_out, flush=True)
          line = '\t'.join([geneid, geneinfo])
          print(line, file=fh_out3, flush=True)

  fh_out.close()
  fh_out2.close()
  fh_out3.close()
