#!/usr/bin/env pytom
import os
import subprocess

# change working directory to tests folder
tests_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(tests_dir)

tests = [test for test in os.listdir('./')
         if test.startswith('test_') and test.endswith('Test.py')]

for test in tests:
    print('\n\nRunning ', test)
    data = subprocess.run(['pytom', test], text=True, capture_output=True).stdout
    if 'FAILED' in data:
        print(f'errors in running {test}') 
