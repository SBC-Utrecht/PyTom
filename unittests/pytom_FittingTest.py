'''
Created on Nov 8, 2010

@author: hrabe
'''

import unittest;

class pytom_FittingTest(unittest.TestCase):
    
    def RefinementFitting_Test(self):
        
        from pytom.fitting.RefinementFitting import refinementFitting
        from pytom.fitting.RefinementFittingStructures import RefinementFittingJob
        from pytom.basic.structures import Particle,Reference,Mask
        from pytom.angles.fromFile import AngleListFromEM
        from pytom.tools.files import getPytomPath
        from pytom.angles.quaternions import Quaternion
        import os
        from pytom_volume import vol
        
        v = vol(32,32,32)
        v.setAll(0)
        v.setV(16,16,16,1)
        
        v.write('fitting_test_volume.em')
        m = vol(32,32,32)
        m.setAll(1)
        m.write('fitting_test_mask.em')
        p = Particle('./fitting_test_volume.em')
        r = Reference('./fitting_test_volume.em')
        m = Mask('./fitting_test_mask.em')
        
        a = AngleListFromEM(getPytomPath() + '/angles/angleLists/angles_000.em',refinementParameters=[0,10])   
        
        job = RefinementFittingJob(p,r,m,a,0,2)

        peak = refinementFitting(job,10,False)
        
        os.system('rm ./fitting_test_volume.em')
        os.system('rm ./fitting_test_mask.em')
        os.system('rm ./referenceOriented.em')