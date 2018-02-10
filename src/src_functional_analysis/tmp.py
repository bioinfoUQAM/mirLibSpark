def binary ():
  '''
  127*.75=95.25
  '''
  infile = 'tmp.txt'
  with open (infile, 'r') as fh:
    for i in fh:
      i = int(i.rstrip('\n'))
      if i > 95: print(75)
      elif i == 1: print(1)
      elif i == 0: print(0)
      else: print(2)
	  
def seqconfidence ():
  infile = 'tmp.txt'
  with open (infile, 'r') as fh:
    for i in fh:
      confid = i.split('_')[0]
      if confid == 'miRB': print('H')
      elif confid == 'low': print('L')
      elif confid == 'novel': print('N')

def trimID ():
  infile = 'tmp.txt'
  with open (infile, 'r') as fh:
    for i in fh:
      id = i.rstrip('\n').split('_')[1]
      print(id)
  
  
  
#binary()
#seqconfidence ()
trimID ()
      
      