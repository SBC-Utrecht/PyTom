"""
test auto-focused classification
"""
import unittest
from pytom.tompy.mpi import MPI


class pytom_MyFunctionTest(unittest.TestCase):

    def setUp(self):
        from helper_functions import create_RandomParticleList

        self.reffile = './testData/ribo.em'
        self.pl_filename = 'pl.xml'
        self.pdir = './testparticles' 
        self.pl = create_RandomParticleList( reffile=self.reffile, pl_filename=self.pl_filename, 
                  pdir=self.pdir, nparticles=10)

        # set parameters for ACL
        self.settings = {}
        self.settings["ncluster"] = 2
        self.settings["frequency"] = 10
        self.settings["binning"] = 4
        self.settings["niteration"] = 2
        #self.settings["fixed_frequency"] = True
        #self.settings["offset"] = None
        #self.settings["mask"] = options.mask
        #self.settings["fmask"] = options.fmask
        #self.settings["mask"] = None
        #self.settings["fmask"] = None
        #self.settings["dispersion"] = None
        #self.settings["external"] = None
        #self.settings["resume"] = options.resume
        #self.settings["resume"] = None
        #self.settings["sigma"] = None
        #self.settings["threshold"] = 0.4
        #self.settings["noise"] = None
        #self.settings["noalign"] = options.noalign
        #self.settings["noalign"] = None
        #self.settings['output_directory'] = './'
        #self.settings["noise"] = 0.1

    #def tearDown(self):
        #self.cleanUp()

    def cleanUp(self):
        from helper_functions import cleanUp_RandomParticleList

        #print('cleaning up AutoFocus files')
        for iclass in range(0, self.settings["ncluster"]):
            tline = 'initial_'+str(iclass)+'.em'
            self.remove_file( filename=tline)
        for ii in range(0,self.settings["niteration"]):
            tline=('classified_pl_iter'+str(ii)+'.xml')
            self.remove_file( filename=tline)
            for iclass in range(0, self.settings["ncluster"]):
                tline=('iter'+str(ii)+'_class'+str(iclass)+'.em')
                self.remove_file( filename=tline)
                tline=('iter'+str(ii)+'_class'+str(iclass)+'_wedge.em')
                self.remove_file( filename=tline)
                for jclass in range(iclass+1, self.settings["ncluster"]):
                    tline=('iter'+str(ii)+'_dmap_'+str(iclass)+'_'+str(jclass)+'.em')
                    self.remove_file( filename=tline)
        cleanUp_RandomParticleList( pl_filename=self.pl_filename, pdir=self.pdir)

    def remove_file(self, filename):
        from os import remove
        try:
            remove(filename)
        except:
            print('file ',filename,' already gone')
                        
    def test_MissingParameters(self):
        """
        call AC without required parameters - check that error is raised
        """
        import os

        cmd = 'mpirun -np 2 ../bin/pytom ../classification/auto_focus_classify.py'
        cmd = cmd + ' -p ' + self.pl_filename
        #cmd = cmd + ' -k ' + str(self.settings["ncluster"])
        #cmd = cmd + ' -f ' + str(self.settings["frequency"])
        try:
            e = os.system(cmd)
        except:
            print('hello')
        self.assertTrue( e == 256, 'no error raised!')


    def test_Binning(self):
        """
        test implementation of binning functionality
        """
        import os

        cmd = 'mpirun -np 2 ../bin/pytom ../classification/auto_focus_classify.py'
        cmd = cmd + ' -p ' + self.pl_filename
        cmd = cmd + ' -k ' + str(self.settings["ncluster"])
        cmd = cmd + ' -f ' + str(self.settings["frequency"])
        cmd = cmd + ' -b ' + str(self.settings["binning"])
        cmd = cmd + ' -i ' + str(self.settings["niteration"])
        print(cmd)
        os.system(cmd)
        self.cleanUp()

        #cmd = cmd + '-o ' + self.settings["ncluster"]
        #cmd = cmd + '-a ' 

        #self.assertAlmostEqual(first=myResult, second=1., places=3, 
        #     msg='the result is not what it is supposed to be')
        #self.assertTrue( result == something, 'wrong result')

        
if __name__ == '__main__':
    #mpi = MPI()
    #mpi.begin()
    unittest.main()
