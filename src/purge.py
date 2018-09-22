import os

cmd = 'rm -fr ../output/*'
os.system(cmd)

cmd = 'rm -fr ../tmp/*'
os.system(cmd)


#cmd = 'su';os.system(cmd)
#cmd = 'cloudera';os.system(cmd)
#cmd = 'free -m';os.system(cmd)
#cmd = 'sync';os.system(cmd)
#cmd = 'echo 1 > /proc/sys/vm/drop_caches';os.system(cmd)
#cmd = 'echo 2 > /proc/sys/vm/drop_caches';os.system(cmd)
#cmd = 'echo 3 > /proc/sys/vm/drop_caches';os.system(cmd)
#cmd = 'free -m';os.system(cmd)
#cmd = 'exit';os.system(cmd)


cmd = 'free -m';os.system(cmd)
cmd = 'sudo sysctl -w vm.drop_caches=3 vm.drop_caches=0';os.system(cmd)
cmd = 'sudo dd if=/dev/sda bs=1M of=/dev/null count=1k';os.system(cmd)
cmd = 'sudo dd if=/dev/sdb bs=1M of=/dev/null count=1k';os.system(cmd)
cmd = 'sudo dd if=/dev/sdc bs=1M of=/dev/null count=1k';os.system(cmd)
cmd = 'sudo dd if=/dev/sdd bs=1M of=/dev/null count=1k';os.system(cmd)
cmd = 'free -m';os.system(cmd)

