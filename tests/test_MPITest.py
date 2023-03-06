"""
Created on May 24, 2013

@author: thrabe
"""
import unittest
import os


class pytom_MPITest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from pytom.tools.files import dump2TextFile

        script_pytom_mpi = 'from pytom.lib.pytom_mpi import rank,init,isInitialised,finalise \n'
        script_pytom_mpi += 'import os\n'
        script_pytom_mpi += 'init()\n'
        script_pytom_mpi += 'os.system("echo " +str(rank())+ " > " +str(rank())+ ".rnkPyTom" )\n'
        script_pytom_mpi += 'finalise()\n'

        dump2TextFile('script_pytom_mpi.py', script_pytom_mpi)

        script_pytom_mpi4py = 'from pytom.agnostic.mpi import MPI\n'
        script_pytom_mpi4py += 'import os\n'
        script_pytom_mpi4py += 'mpi = MPI()\n'
        script_pytom_mpi4py += 'os.system("echo " + str(mpi.rank) + " > " + str(mpi.rank) + ".rnkPyTom")\n'
        script_pytom_mpi4py += 'mpi.begin()\n'
        # this will put the follower in standby and only allow the leader to execute
        script_pytom_mpi4py += 'os.system("echo " + str(mpi.rank) + " > leader" + str(mpi.rank) + str(mpi.size) + '
        script_pytom_mpi4py += '".mpiLeader")\n'
        script_pytom_mpi4py += 'mpi.end()\n'

        dump2TextFile('script_pytom_mpi4py.py', script_pytom_mpi4py)

    @classmethod
    def tearDownClass(cls):
        os.system('rm -f script_pytom_mpi.py')
        os.system('rm -f script_pytom_mpi4py.py')

    def tearDown(self):
        os.system('rm -f *.rnkPyTom')  # remove also in case any one of the tests failed
        os.system('rm -f *.mpiLeader')
        
    def test_mpi(self):
        os.system('mpirun -n 2 python script_pytom_mpi.py')
        ranks = sorted([int(f.strip('.rnkPyTom')) for f in os.listdir('./') if f.endswith('.rnkPyTom')])
        self.assertTrue(len(ranks) == 2 and ranks[0] < ranks[1], msg='mpi ranks error in pytom.lib.pytom_mpi')

    def test_mpi4py(self):
        os.system('mpirun -n 2 python script_pytom_mpi4py.py')
        ranks = sorted([int(f.strip('.rnkPyTom')) for f in os.listdir('./') if f.endswith('.rnkPyTom')])
        self.assertTrue(len(ranks) == 2 and ranks[0] < ranks[1],
                        msg='mpi ranks error in pytom.agnostic.mpi via mpi4py')
        leader_msg = [f for f in os.listdir('./') if f.endswith('.mpiLeader')]
        self.assertTrue(len(leader_msg) == 1 and leader_msg[0] == 'leader02.mpiLeader',
                        msg='worker was not correctly put on standby by leader')


if __name__ == "__main__":
    unittest.main()
