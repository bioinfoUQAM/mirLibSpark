import time, os

for i in range(1, 81):
  cmd = 'sbatch ' + str(i) + '.sh';print(cmd);os.system(cmd)
  time.sleep (5)
