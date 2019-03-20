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
  '''
  When the positive and negative classes are unbalanced, F1 can not be used. Instead, use MCC.
  
  The MCC value is between -1 and 1.
  0 : random
  -1: disagree
  1 : agree
  '''
  up = tp*tn - fp*fn
  dn = sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  return round (up/dn, 3)
  
def report (TP, TN, FP, FN):
  f1 = calculate_F1_score (TP, TN, FP, FN)
  acc = calculate_accuracy (TP, TN, FP, FN)
  sens = sensitivity (TP, FN)
  mcc = Matthews_correlation_coefficient (TP, TN, FP, FN)

  print('TP\tTN\tFP\tFN\tF1\tAcc\tMCC\tsensitivity')
  data = [str(x) for x in [TP, TN, FP, FN, f1, acc, mcc, sens]]
  line = '\t'.join(data)
  print(line)
  print()

msg = 'mirlibspark';print(msg)
#== simulation ==
TP = 22
FP = 0
msg = 'simulation';print(msg)
report (TP, 1000 - FP, FP, 100 - TP)

#== 100.txt ==
TP = 29
FP = 343
negset = 1838879 - 141
msg = 'case study';print(msg)
report ( TP, (negset-FP), FP, (141-TP) )
#report (TP, TN, FP, FN)



