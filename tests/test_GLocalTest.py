"""
test GLocal Alignment
"""
import unittest

class pytom_GLocalTest(unittest.TestCase):

    def setUp(self):
        from helper_functions import create_RandomParticleList, installdir
        from pytom.tompy.io import read_size
        from pytom_volume import vol, initSphere

        self.installdir = installdir
        self.reffile = f'./testData/ribo.em'
        self.pl_filename = 'pl.xml'
        self.pdir = f'./testparticles'
        self.pl = create_RandomParticleList( reffile=self.reffile, pl_filename=self.pl_filename, 
                  pdir=self.pdir, nparticles=10)

        # set parameters for GLocal
        self.settings = {}
        self.settings["binning"] = 4
        self.settings["niteration"] = 1
        self.settings["mask"] = f'./testData/ribo_mask.em'
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

        for ii in range(0, self.settings["niteration"]):
            fname = str(ii)+'-All.mrc'
            self.remove_file( filename=fname)
            fname = str(ii)+'-AllFiltered*.mrc'
            system('rm '+fname)
            fname = str(ii)+'-Even.mrc'
            self.remove_file( filename=fname)
            fname = str(ii)+'-EvenFiltered.mrc'
            self.remove_file( filename=fname)
            fname = str(ii)+'-EvenFiltered-PreWedge.mrc'
            self.remove_file( filename=fname)
            fname = str(ii)+'-EvenFiltered-WedgeSumUnscaled.mrc'
            self.remove_file( filename=fname)
            fname = str(ii)+'-Odd.mrc'
            self.remove_file( filename=fname)
            fname = str(ii)+'-OddFiltered.mrc'
            self.remove_file( filename=fname)
            fname = str(ii)+'-OddFiltered-PreWedge.mrc'
            self.remove_file( filename=fname)
            fname = str(ii)+'-OddFiltered-WedgeSumUnscaled.mrc'
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
        self.remove_file( filename='average-Final.mrc')
        self.remove_file( filename='average-Final-PreWedge.mrc')
        self.remove_file( filename='average-Final-WedgeSumUnscaled.mrc')
        self.remove_file( filename='FSC-Final.dat')
        self.remove_file( filename='myGLocal.xml')
        self.remove_file( filename='average-Final-Even.mrc')
        self.remove_file( filename='average-Final-Even-PreWedge.mrc')
        self.remove_file( filename='average-Final-Even-WedgeSumUnscaled.mrc')
        self.remove_file( filename='average-Final-Odd.mrc')
        self.remove_file( filename='average-Final-Odd-PreWedge.mrc')
        self.remove_file( filename='average-Final-Odd-WedgeSumUnscaled.mrc')
        fname='average-FinalFiltered_*.mrc'
        system('rm '+fname)
        system(f'rm {self.pl_filename}')
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

        cmd = f'mpirun -np 2 {self.installdir}/bin/pytom {self.installdir}/bin/GLocalJob.py'
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

        cmd = f'mpirun -np 2 {self.installdir}/bin/pytom {self.installdir}/bin/GLocalJob.py'
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
