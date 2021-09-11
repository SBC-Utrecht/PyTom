#!/usr/bin/env pytom
import os

tests = [test for test in os.listdir('./') if test.startswith('test_') and test.endswith('Test.py')]

for test in tests:
    print('\n\nRunning ', test)
    data = os.popen(f'pytom {test}').read()
    if 'FAILED' in data:
        print(f'errors in running {test}') 
