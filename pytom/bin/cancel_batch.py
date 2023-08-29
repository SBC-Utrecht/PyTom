#!/usr/bin/env pytom

import os, sys
import subprocess
import getpass
from shutil import which

if which('squeue') is None or which('scancel') is None:
    print('SLURM commands not available in PATH, exiting')
    exit()

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

username = getpass.getuser()
print(username, start, end)
running = [int(line) for line in 
        subprocess.run(["squeue -u " + username + " | awk '{print $1}'"], 
            shell=True, capture_output=True, text=True
            ).stdout.splitlines()[1:] if line]

jobs = list(range(start, end+1))
if jobs and running:
    # Make sure that both of these lists are not empty
    start = max(min(jobs), min(running))

    for i in range(start, end+1):
        print(i)
        if i in running: os.system('scancel {}'.format(i))
