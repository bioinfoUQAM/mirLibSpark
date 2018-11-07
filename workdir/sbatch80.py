import time, os

for i in range(33, 35):
  cmd = 'sbatch ' + str(i) + '.sh';print(cmd);os.system(cmd)
  time.sleep (5)
