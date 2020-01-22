"""
test GLocal Alignment
"""
import unittest

class pytom_GLocalTest(unittest.TestCase):

    def setUp(self):
        from helper_functions import create_RandomParticleList
        from pytom.tompy.io import read_size
        from pytom_volume import vol, initSphere

        self.reffile = './testData/ribo.em'
        self.pl_filename = 'pl.xml'
        self.pdir = './testparticles' 
        #self.pl = create_RandomParticleList( reffile=self.reffile, pl_filename=self.pl_filename, 
        #          pdir=self.pdir, nparticles=10)

        # set parameters for GLocal
        self.settings = {}
        self.settings["binning"] = 4
        self.settings["niteration"] = 1
        self.settings["mask"] = './testData/ribo_mask.em'
        dims = read_size(self.reffile)
        maskvol = vol(int(dims[0]), int(dims[1]), int(dims[2]))
        initSphere(maskvol, 30,5, 0, int(dims[0]/2), int(dims[1]/2), int(dims[2]/2))
        maskvol.write(self.settings["mask"])
        self.settings["destination"] = './'
        #self.settings["score"] = 'nxcf'
        self.settings["score"] = 'flcf'
        self.settings["pixelsize"] = 2.
        self.settings["diameter"] = 250.
        self.settings["job"] = './myGLocal.xml'

    #def tearDown(self):
        #self.cleanUp()

    def cleanUp(self):
        """
        check that files are written and remove them
        """
        from helper_functions import cleanUp_RandomParticleList
        from os import system

        for ii in range(0, self.settings["niteration"]+1):
            fname = str(ii)+'-All.em'
            self.remove_file( filename=fname)
            fname = str(ii)+'-AllFiltered*.em'
            system('rm '+fname)
            fname = str(ii)+'-Even.em'
            self.remove_file( filename=fname)
            fname = str(ii)+'-EvenFiltered.em'
            self.remove_file( filename=fname)
            fname = str(ii)+'-EvenFiltered-PreWedge.em'
            self.remove_file( filename=fname)
            fname = str(ii)+'-EvenFiltered-WedgeSumUnscaled.em'
            self.remove_file( filename=fname)
            fname = str(ii)+'-Odd.em'
            self.remove_file( filename=fname)
            fname = str(ii)+'-OddFiltered.em'
            self.remove_file( filename=fname)
            fname = str(ii)+'-OddFiltered-PreWedge.em'
            self.remove_file( filename=fname)
            fname = str(ii)+'-OddFiltered-WedgeSumUnscaled.em'
            self.remove_file( filename=fname)
            fname = str(ii)+'-FSC.dat'
            self.remove_file( filename=fname)
            fname = str(ii)+'-Filter.dat'
            self.remove_file( filename=fname)
            fname = str(ii)+'-GLocalAlignmentJob.xml'
            self.remove_file( filename=fname)
            fname = str(ii)+'-ParticleListEven.xml'
            self.remove_file( filename=fname)
            fname = str(ii)+'-ParticleListOdd.xml'
            self.remove_file( filename=fname)
            fname = str(ii)+'-ParticleList.xml'
            self.remove_file( filename=fname)
        self.remove_file( filename='average-Final.em')
        self.remove_file( filename='average-Final-PreWedge.em')
        self.remove_file( filename='average-Final-WedgeSumUnscaled.em')
        self.remove_file( filename='FSC-Final.dat')
        self.remove_file( filename='myGLocal.xml')
        self.remove_file( filename='average-Final-Even.em')
        self.remove_file( filename='average-Final-Even-PreWedge.em')
        self.remove_file( filename='average-Final-Even-WedgeSumUnscaled.em')
        self.remove_file( filename='average-Final-Odd.em')
        self.remove_file( filename='average-Final-Odd-PreWedge.em')
        self.remove_file( filename='average-Final-Odd-WedgeSumUnscaled.em')
        fname='average-FinalFiltered_*.em'
        system('rm '+fname)

        #cleanUp_RandomParticleList( pl_filename=self.pl_filename, pdir=self.pdir)

    def remove_file(self, filename):
        """
        assert that file exists end remove it
        """
        from os import remove
        from os import path

        filecheck = path.exists(filename)
        if not filecheck:
            print(filename+" does not exist")
        self.assertTrue( filecheck, msg="file "+filename+" does not exist")
        if filecheck:
            remove(filename)


    def xtest_Score(self):
        """
        test implementation of different scores
        """
        import os

        cmd = 'mpirun -np 2 ../bin/pytom ../bin/GLocalJob.py'
        cmd = cmd + ' -p ' + self.pl_filename
        cmd = cmd + ' -m ' + str(self.settings["mask"])
        cmd = cmd + ' -b ' + str(self.settings["binning"])
        cmd = cmd + ' -n ' + str(self.settings["niteration"])
        cmd = cmd + ' -d ' + str(self.settings["destination"])
        cmd = cmd + ' -s ' + str(self.settings["score"])
        cmd = cmd + ' --pixelSize ' + str(self.settings["pixelsize"])
        cmd = cmd + ' --particleDiameter ' + str(self.settings["diameter"])
        cmd = cmd + ' -j ' + str(self.settings["job"])
        #cmd = cmd + ' -r ' + str(self.settings["reference"])
        print(cmd)
        os.system(cmd)
        self.cleanUp()

    def test_nXcf(self):
        """
        test implementation of different scores
        """
        import os

        self.settings["score"] = 'nxcf'

        cmd = 'mpirun -np 2 ../bin/pytom ../bin/GLocalJob.py'
        cmd = cmd + ' -p ' + self.pl_filename
        cmd = cmd + ' -m ' + str(self.settings["mask"])
        cmd = cmd + ' -b ' + str(self.settings["binning"])
        cmd = cmd + ' -n ' + str(self.settings["niteration"])
        cmd = cmd + ' -d ' + str(self.settings["destination"])
        cmd = cmd + ' -s ' + str(self.settings["score"])
        cmd = cmd + ' --pixelSize ' + str(self.settings["pixelsize"])
        cmd = cmd + ' --particleDiameter ' + str(self.settings["diameter"])
        cmd = cmd + ' -j ' + str(self.settings["job"])
        #cmd = cmd + ' -r ' + str(self.settings["reference"])
        print(cmd)
        os.system(cmd)
        self.cleanUp()

        #self.assertAlmostEqual(first=myResult, second=1., places=3, 
        #     msg='the result is not what it is supposed to be')
        #self.assertTrue( result == something, 'wrong result')


        
if __name__ == '__main__':
    unittest.main()
