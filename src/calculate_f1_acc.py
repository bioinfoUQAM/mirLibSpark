def calculate_F1_score (tp, tn, fp, fn):
  #= p: precision; r: recall
  p = tp / (tp + fp)
  r = tp / (tp + fn)
  f1 = 2 * p * r / (p + r)
  return round (f1, 3)

def calculate_accuracy (tp, tn, fp, fn):
  acc = (tp + tn)/(tp + tn + fp + fn)
  return round (acc, 3)


TP = 73
TN = 999
FP = 1
FN = 27

  
f1 = calculate_F1_score (TP, TN, FP, FN)
acc = calculate_accuracy (TP, TN, FP, FN)

print(f1)
print(acc)

