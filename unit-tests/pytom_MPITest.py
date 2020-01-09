#!/usr/bin/env pytom
'''
Created on May 24, 2013

@author: thrabe
'''
import unittest


class pytom_MPITest(unittest.TestCase):

    def is_exe(self, file):
        import os
        return os.path.isfile(file) and os.access(file, os.X_OK)

    def setUp(self):
        from pytom.tools.files import dump2TextFile
        #from pytom_mpi import rank
        #exeNames = ['openmpirun','mpirun']
        
        #runMPI  = [exeNames[i] for i in range(len(exeNames)) if self.is_exe(exeNames[i])]
        runMPI = 'mpirun'
        print(runMPI)
        self.testCommand = runMPI + ' -c 2 pytom testScript.py'
         
        testScript  = 'from pytom_mpi import rank,init,isInitialized,finalise \n'
        testScript += 'import os\n'
        testScript += 'init()\n'
        testScript += 'os.system("echo " +str(rank())+ " > " +str(rank())+ ".rnkPyTom" )\n'
        testScript += 'finalise()\n'

        dump2TextFile('testScript.py',testScript)


    def tearDown(self):
        import os
        
        os.system('rm -f *rnkPyTom')
        os.system('rm -f testScript.py')
        
    def test_mpi(self):
        import os
        os.system(self.testCommand)

        pass

    

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
