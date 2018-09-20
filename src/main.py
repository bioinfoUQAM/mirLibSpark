import os

#cmd = 'python purge.py'
#os.system (cmd)

cmd = 'git pull origin smooth'
os.system (cmd)

cmd = 'time spark-submit mirLibPipeline.py ../paramfile.txt 2>/dev/null'
os.system(cmd)
