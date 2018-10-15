'''
Amine M Remita
Chao-Jung Wu

2018-10-05
'''
import math
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests


def __calculate_z_score (N1, N2, n1, n2):
  '''  # Kal et al 1999  '''
  p0 = (n1 + n2)/(N1 + N2)
  p1 = n1/N1
  p2 = n2/N2
  denom = math.sqrt(p0 * (1-p0) * ( (1/N1) + (1/N2)))
  if denom == 0: zscore = 0
  else: zscore = (p1-p2)/denom
  return zscore

def __calculate_population_sum (data):
  N1, N2 = 0, 0
  for i in data:   
    a = float(i[1])
    b = float(i[2])
    N1 += a
    N2 += b
  return N1, N2

def __calculate_BH_p_values (DATA):
  P_VALUES = []
  for n in range(len(DATA)):
    i = DATA[n]
    p_value = i[5]
    P_VALUES.append(p_value)

  #= https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html
  #= reject (array, boolean) – true for hypothesis that can be rejected for given alpha
  #= pvals_corrected (array) – p-values corrected for multiple tests
  #= alphacSidak (float) – corrected alpha for Sidak method
  #= alphacBonf (float) – corrected alpha for Bonferroni method
  reject, pvals_corrected, alphacSidak, alphacBonf = multipletests(P_VALUES)
  return pvals_corrected




#infile = '../input/demo_diff_input.txt'
#with open (infile, 'r') as fh: data = [ x.rstrip('\n').split('\t') for x in fh.readlines()]


#title = 'Sequence	Norstar_y2013_5	Norstar_y2013_1	Fold_change	Z_score	p_value	BH_p_value	Diff_exp'.split('\t')

def main (data, title):
  N1, N2 = __calculate_population_sum (data[1:])
  DATA = []
  DATA.append(title)
  for n in range(1, len(data)):
    i = data[n]
    seq = i[0]    
    a = float(i[1])
    b = float(i[2])

    foldchange = ( a + 1 ) / ( b + 1 )
    zscore = __calculate_z_score(N1, N2, a, b)
    p_value = norm.cdf(zscore) #= https://docs.scipy.org/doc/scipy-0.16.1/reference/generated/scipy.stats.norm.html
    i.append(foldchange)
    i.append(zscore)
    i.append(p_value)
    DATA.append(i)

  list_BH_p_values = __calculate_BH_p_values (DATA[1:])

  for n in range(1, len(data)):
    DATA[n].append(list_BH_p_values[n-1])
  return DATA
  
  #for i in DATA: print( '\t'.join([str(x) for x in i]) )
  



  
  
  
  
  
