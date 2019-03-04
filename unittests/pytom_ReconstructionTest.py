'''
Created on Mar 2, 2011

@author: luiskuhn
'''
import unittest;

class pytom_ReconstructionTest(unittest.TestCase):
    
    def weightedBackprojection_Test(self):
        
        from pytom_volume import vol,backProject,mean
        
        from pytom.basic.normalise import mean0std1
        from pytom.reconstruction.weightedBackprojection import getWeightedProjectionCube

        from pytom.tools.functions import frange
        import math

        vol_src = vol(100,100,100)
        vol_src.setV(1,50,50,50)

        vol_src_dim_x = vol_src.sizeX()
        vol_src_dim_y = vol_src.sizeY()
        vol_src_dim_z = vol_src.sizeZ()

        angle_range = frange(-90, 90, 0.5)
        angleNo = len(angle_range)



        vol_bp = vol(vol_src.sizeX(), vol_src.sizeY(), vol_src.sizeZ())
        vol_bp.setAll(0.0)

        vol_phi = vol(1, 1, len(angle_range))
        vol_phi.setAll(0.0)

        vol_the = vol(1, 1, len(angle_range))
        vol_the.setAll(0.0)

        vol_offset = vol(3, len(self) ,1)
        vol_offset.setAll(0.0)

        vol_offsetProjections = vol(1, 2, len(self))
        vol_offsetProjections.setAll(0.0)

        for i in range(angleNo):
            angle = angle_range[i]    
            vol_the.setV(-angle, 0, 0, i)

        vol_img = getWeightedProjectionCube(vol_src, angle_range)


        backProject(vol_img, vol_bp, vol_phi, vol_the, vol_offset,vol_offsetProjections)
        mean0std1(vol_bp)

        
        eps = 5

        assert math.fabs(mean(subvolume(vol_bp, 40, 40, 40, 20, 20, 20))-vol_src.getV(50, 50, 50)) < eps
        
    def goldReconstruction_Test(self):
        
        from pytom_volume import read,mean,variance
        from pytom.tools.files import getPytomPath
        from pytom.reconstruction.reconstructionStructures import ProjectionList
        from pytom.basic.structures import Particle
        import os
        cwd = os.getcwd()
        
        p = Particle('./testData/reconstruction/p1_new.em')

        projList = ProjectionList()
        projList.fromXMLFile('./testData/reconstruction/p1/projectionList.xml')
        projList.reconstructVolumes( particles=[p], cubeSize=100, binning=1, applyWeighting=False,
                                     showProgressBar=False, verbose=False, preScale=1, postScale=1)
        
        v1 = read('./testData/reconstruction/p1.em')
        v2 = read('./testData/reconstruction/p1_new.em')
        
        assert ( abs(mean(v1) - mean(v2)) < 0.0000000001) and (abs(variance(v1, True) - variance(v2, True)) < 0.0000000001)
        
        os.remove('./testData/reconstruction/p1_new.em')
        os.chdir(cwd)
    
        
