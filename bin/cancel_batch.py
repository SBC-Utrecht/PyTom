#!/usr/bin/env pytom

import os, sys
import getpass

if len(sys.argv[1:]) == 2:
    try:
       start,end = map(int,sys.argv[1:])
    except:
       raise Exception('invalid input parameters')
elif len(sys.argv[1:]) == 1:
    try:
       start = int(sys.argv[1])
       end = start+1
    except:
       raise Exception('invalid input parameters')
else:
    raise Exception('invalid input parameters')

#username = os.popen('whoami').read()[:-1]
username = getpass.getuser()
print(username, start, end)
running = [int(line) for line in os.popen("squeue -u " + username + " | awk '{print $1}'").readlines()[1:] if line]

jobs = list(range(start, end+1))
start = max(min(jobs), min(running))

for i in range(start, end+1):
    print(i)
    if i in running: os.system('scancel {}'.format(i))
