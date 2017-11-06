def calculate_F1_score (tp, tn, fp, fn):
  #= p: precision; r: recall
  p = tp / float(tp + fp)
  r = tp / float(tp + fn)
  f1 = 2 * p * r / (p + r)
  return round (f1, 3)

def calculate_accuracy (tp, tn, fp, fn):
  acc = (tp + tn)/ float(tp + tn + fp + fn)
  return round (acc, 3)


TP = 82
TN = 990
FP = 10
FN = 18

  
f1 = calculate_F1_score (TP, TN, FP, FN)
acc = calculate_accuracy (TP, TN, FP, FN)

print(f1)
print(acc)

