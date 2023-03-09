"""
test auto-focused classification
"""
import unittest
import os

class pytom_MyFunctionTest(unittest.TestCase):

    def setUp(self):
        from helper_functions import create_RandomParticleList, installdir

        self.installdir = installdir
        self.reffile = f'../testData/ribo.em'
        self.pdir = f'./testparticles'


        # set parameters for ACL
        self.settings = {}
        self.settings["frequency"] = 40
        self.settings["outputDirectory"] = 'outputCCC'
        if not os.path.exists(self.settings["outputDirectory"]):
            os.mkdir(self.settings["outputDirectory"])

        self.pl_filename = os.path.join(self.settings['outputDirectory'], 'pl.xml')

        self.pl = create_RandomParticleList(reffile=self.reffile, pl_filename=self.pl_filename,
                                            pdir=self.pdir, nparticles=10)
        #self.settings["fixed_frequency"] = True
        #self.settings["offset"] = None
        #self.settings["mask"] = options.mask

        self.settings["mask"] = f'../testData/mask_45.em'


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
        """
        check that files are written and remove them
        """
        from helper_functions import cleanUp_RandomParticleList
        import os
        # for iclass in range(0, self.settings["ncluster"]):
        #     tline = 'initial_'+str(iclass)+'.em'
        #     self.remove_file( filename=tline)
        # for ii in range(0,self.settings["niteration"]):
        #     tline=('classified_pl_iter'+str(ii)+'.xml')
        #     self.remove_file( filename=tline)
        #     for iclass in range(0, self.settings["ncluster"]):
        #         tline=('iter'+str(ii)+'_class'+str(iclass)+'.em')
        #         self.remove_file( filename=tline)
        #         tline=('iter'+str(ii)+'_class'+str(iclass)+'_wedge.em')
        #         self.remove_file( filename=tline)
        #         for jclass in range(iclass+1, self.settings["ncluster"]):
        #             tline=('iter'+str(ii)+'_dmap_'+str(jclass)+'_'+str(iclass)+'.em')
        #             self.remove_file( filename=tline)
        #             tline=('iter'+str(ii)+'_dmap_'+str(iclass)+'_'+str(jclass)+'.em')
        #             self.remove_file( filename=tline)

        self.remove_file( filename=f'{self.settings["outputDirectory"]}/correlation_matrix.csv')

        cleanUp_RandomParticleList( pl_filename=self.pl_filename, pdir=self.pdir)
        os.removedirs(self.settings["outputDirectory"])

    def remove_file(self, filename):
        """
        assert that file exists end remove it
        """
        from os import remove
        from os import path

        filecheck = path.exists(filename)
        self.assertTrue( filecheck, msg="file "+filename+" does not exist")
        if filecheck:
            remove(filename)
        #else:
        #    print('file ',filename,' already gone')


    def test_CCC_CPU(self):
        """
        test implementation of binning functionality
        """
        import os

        cmd = f'mpirun -np 2 pytom {self.installdir}/classification/calculate_correlation_matrix.py'
        cmd = cmd + ' -p ' + self.pl_filename
        cmd = cmd + ' -f ' + str(self.settings["frequency"])
        cmd = cmd + ' -m ' + str(self.settings["mask"])
        cmd = cmd + ' -o ' + str(self.settings["outputDirectory"])
        print(cmd)
        os.system(cmd)
        self.cleanUp()

    def CCC_GPU(self):
        """
        test gpu implementation
        """
        import os

        cmd = f'mpirun -np 2 pytom {self.installdir}/classification/calculate_correlation_matrix.py'
        cmd = cmd + ' -p ' + self.pl_filename
        cmd = cmd + ' -f ' + str(self.settings["frequency"])
        cmd = cmd + ' -m ' + str(self.settings["mask"])
        cmd = cmd + ' -o ' + str(self.settings["outputDirectory"])
        cmd += ' -g 0'
        print(cmd)
        os.system(cmd)
        os.system(f'rm -rf {self.settings["outputDirectory"]}')

        self.cleanUp()



        
if __name__ == '__main__':
    unittest.main()
