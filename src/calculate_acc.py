'''
Author: Chao-Jung Wu
Date: 2017-11-07
'''
from __future__ import print_function
from math import sqrt


def calculate_F1_score (tp, tn, fp, fn):
  #= p: precision; r: recall
  p = tp / float(tp + fp)
  r = tp / float(tp + fn)
  f1 = 2 * p * r / (p + r)
  return round (f1, 3)

def calculate_accuracy (tp, tn, fp, fn):
  acc = (tp + tn)/ float(tp + tn + fp + fn)
  return round (acc, 3)

def sensitivity (tp, fn):
  sens = tp / float(tp + fn)
  return round (sens, 3)

def Matthews_correlation_coefficient (tp, tn, fp, fn):
  up = tp*tn - fp*fn
  dn = sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  return round (up/dn, 3)


TP = 82
TN = 990
FP = 10
FN = 18

#===========
TP = 84
FP = 9

FN = 100 - TP
TN = 1000 - FP
#===========

##===========
## 100.txt
#mature = 71
#FP = 52
#TP = mature - FP
#
#FN = 141 - TP
#TN = 1838738 - FP
#===========
  
f1 = calculate_F1_score (TP, TN, FP, FN)
acc = calculate_accuracy (TP, TN, FP, FN)
sens = sensitivity (TP, FN)
mcc = Matthews_correlation_coefficient (TP, TN, FP, FN)


print('TP\tTN\tFP\tFN\tF1\tAccuracy')
data = [str(x) for x in [TP, TN, FP, FN, f1, acc]]
line = '\t'.join(data)
print(line)


print()
print('TP\tTN\tFP\tFN\tF1\tMCC')
data = [str(x) for x in [TP, TN, FP, FN, f1, mcc]]
line = '\t'.join(data)
print(line)
