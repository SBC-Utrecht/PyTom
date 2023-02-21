"""
test GLocal Alignment
"""
import unittest
import os


class pytom_GLocalTest(unittest.TestCase):

    def setUp(self):
        from helper_functions import create_RandomParticleList, installdir
        from pytom.agnostic.io import read_size
        from pytom.lib.pytom_volume import vol, initSphere

        self.installdir = installdir
        self.reffile = f'../testData/ribo.em'
        self.pl_filename = 'pl.xml'
        self.pdir = f'./testparticles'
        self.pl = create_RandomParticleList( reffile=self.reffile, pl_filename=self.pl_filename, 
                  pdir=self.pdir, nparticles=10)

        # set parameters for GLocal
        self.settings = {}
        self.settings["binning"] = 1
        self.settings["niteration"] = 1
        self.settings["mask"] = f'../testData/ribo_mask.em'
        dims = read_size(self.reffile)
        maskvol = vol(int(dims[0]), int(dims[1]), int(dims[2]))
        initSphere(maskvol, 30,5, 0, int(dims[0]/2), int(dims[1]/2), int(dims[2]/2))
        maskvol.write(self.settings["mask"])
        self.settings["destination"] = './'
        self.settings["score"] = 'nxcf'
        self.settings["pixelsize"] = 2.
        self.settings["diameter"] = 250.
        self.settings["job"] = './myGLocal.xml'

    def cleanUp(self):
        """
        check that files are written and remove them
        """

        for ii in range(0, self.settings["niteration"]):
            fname = str(ii)+'-All.em'
            self.remove_file( filename=fname)
            fname = str(ii)+'-AllFiltered*.em'
            os.system('rm '+fname)
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
        os.system('rm '+fname)
        os.system(f'rm {self.pl_filename}')

    def remove_file(self, filename):
        """
        assert that file exists end remove it
        """
        filecheck = os.path.exists(filename)
        if not filecheck:
            print(filename+" does not exist")
        self.assertTrue( filecheck, msg="file "+filename+" does not exist")
        if filecheck:
            os.remove(filename)

    def test_flcf_cpu(self):
        """
        test glocal cpu with fast local correlation function
        """
        cmd = f'mpirun -n 1 GLocalJob.py'
        cmd = cmd + ' -p ' + self.pl_filename
        cmd = cmd + ' -m ' + str(self.settings["mask"])
        cmd = cmd + ' --SphericalMask'
        cmd = cmd + ' -b ' + str(self.settings["binning"])
        cmd = cmd + ' -n ' + str(self.settings["niteration"])
        cmd = cmd + ' -d ' + str(self.settings["destination"])
        cmd = cmd + ' -s ' + 'flcf'
        cmd = cmd + ' --pixelSize ' + str(self.settings["pixelsize"])
        cmd = cmd + ' --particleDiameter ' + str(self.settings["diameter"])
        cmd = cmd + ' -j ' + str(self.settings["job"])
        print(cmd)
        os.system(cmd)
        self.cleanUp()

    def test_flcf_gpu(self):
        """
        test glocal cpu with fast local correlation function
        """
        cmd = f'mpirun -n 1 GLocalJob.py -g 0'
        cmd = cmd + ' -p ' + self.pl_filename
        cmd = cmd + ' -m ' + str(self.settings["mask"])
        cmd = cmd + ' --SphericalMask'
        cmd = cmd + ' -b ' + str(self.settings["binning"])
        cmd = cmd + ' -n ' + str(self.settings["niteration"])
        cmd = cmd + ' -d ' + str(self.settings["destination"])
        cmd = cmd + ' -s ' + 'flcf'
        cmd = cmd + ' --pixelSize ' + str(self.settings["pixelsize"])
        cmd = cmd + ' --particleDiameter ' + str(self.settings["diameter"])
        cmd = cmd + ' -j ' + str(self.settings["job"])
        print(cmd)
        os.system(cmd)
        self.cleanUp()

    def test_nxcf_gpu(self):
        """
        test glocal gpu with normalised cross correlation function (nxcf)
        """
        # had to change number of mpi procs to 1 to comply with running on single gpu
        # BUT: All the docs say that the number of mpi processes should be ( number_of_gpus + 1 )
        cmd = f'mpirun -n 1 GLocalJob.py -g 0'
        cmd = cmd + ' -p ' + self.pl_filename
        cmd = cmd + ' -m ' + str(self.settings["mask"])
        cmd = cmd + ' --SphericalMask'
        cmd = cmd + ' -b ' + str(self.settings["binning"])
        cmd = cmd + ' -n ' + str(self.settings["niteration"])
        cmd = cmd + ' -d ' + str(self.settings["destination"])
        cmd = cmd + ' -s ' + str(self.settings["score"])
        cmd = cmd + ' --pixelSize ' + str(self.settings["pixelsize"])
        cmd = cmd + ' --particleDiameter ' + str(self.settings["diameter"])
        cmd = cmd + ' -j ' + str(self.settings["job"])
        print(cmd)
        os.system(cmd)
        self.cleanUp()

    def test_nxcf_cpu(self):
        """
        test glocal cpu with normalised cross correlation function (nxcf)
        """
        cmd = f'mpirun -n 1 GLocalJob.py '
        cmd = cmd + ' -p ' + self.pl_filename
        cmd = cmd + ' -m ' + str(self.settings["mask"])
        cmd = cmd + ' --SphericalMask'
        cmd = cmd + ' -b ' + str(self.settings["binning"])
        cmd = cmd + ' -n ' + str(self.settings["niteration"])
        cmd = cmd + ' -d ' + str(self.settings["destination"])
        cmd = cmd + ' -s ' + str(self.settings["score"])
        cmd = cmd + ' --pixelSize ' + str(self.settings["pixelsize"])
        cmd = cmd + ' --particleDiameter ' + str(self.settings["diameter"])
        cmd = cmd + ' -j ' + str(self.settings["job"])
        print(cmd)
        os.system(cmd)
        self.cleanUp()

        
if __name__ == '__main__':
    unittest.main()
