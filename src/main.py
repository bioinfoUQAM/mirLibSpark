import os

#cmd = 'python purge.py'
#os.system (cmd)

cmd = 'git pull origin smooth'
os.system (cmd)

#cmd = 'rm -fr ../input/fake_a3.txt'
#os.system(cmd)


cmd = 'cp ../input_samples/fake_a3.txt ../input'
os.system(cmd)


cmd = 'time spark-submit mirLibPipeline.py ../paramfile.txt 2>/dev/null'
os.system(cmd)
