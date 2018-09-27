from __future__ import print_function

infile = '../input_samples/miranda_output_demo.txt'
outfile = '../output/mirna_vs_targetgenes.txt' 

def process_targets(infile, outfile):
  fh_out = open (outfile, 'w')
  print( '\t'.join('mirna, top1TG, top5scoredTG'.split(', ')), file=fh_out, flush=True )
  with open (infile, 'r') as fh: DATA = [x.rstrip('\n') for x in fh.readlines()]

  for i in DATA:
    i=i.split(' [[')
    mirnaseq = i[0]
    i=i[1].split(']] ')
    mirna_uindex = i[1]
    targets=i[0].split('], [')
    targetcollect = []
    score = 1000000000000
    count = 0
    for t in targets:
      target_miranda = t.strip('\'').split('\', \'')
      target = target_miranda[0].split('.')[0]
      score_cur = int(target_miranda[1].split('.')[0])
      if score_cur < score:
        count += 1
        score = score_cur
      if count < 6 and target not in targetcollect:
        targetcollect.append( target + ' ('+ str(score_cur) +')' )
    data = [mirnaseq, targetcollect[0], ','.join(targetcollect)]
    line = '\t'.join(data)
    print(line, file=fh_out, flush=True)

#process_targets(infile, outfile)
