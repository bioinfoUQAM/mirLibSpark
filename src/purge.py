import os

cmd = 'rm -fr ../output/*'
os.system(cmd)

cmd = 'rm -fr ../tmp/*'
os.system(cmd)



cmd = 'free -m'
os.system(cmd)
cmd = 'sync'
os.system(cmd)
cmd = 'echo 1 > /proc/sys/vm/drop_caches'
os.system(cmd)
cmd = 'echo 2 > /proc/sys/vm/drop_caches'
os.system(cmd)
cmd = 'echo 3 > /proc/sys/vm/drop_caches'
os.system(cmd)
cmd = 'free -m'
os.system(cmd)



