'''
Created on Oct 6, 2014

@author: thrabe
'''

import unittest

class pytom_SymmetryTest(unittest.TestCase):
    
    
    
    def HelicalSymmetry_ApplyToParticleTest(self):
        
        from pytom.basic.structures import HelicalSymmetry
        from pytom.basic.files import read
        
        v = read('./testData/helicalSymmetryReference.em')
        vResult = read('./testData/helSymLeftResult.em')
        
        helSymObject = HelicalSymmetry(nfold = 4, symmetryAngle = 360/13.0, repeatShift = 21,isRightSymmetry = False)
        
        v2 = helSymObject.applyToPatricle(v)
        
        assert vResult.equalsTo(v2)
        
    def HelicalSymmetry_ApplyToParticleList(self):
        
        from pytom.basic.structures import HelicalSymmetry,Particle,ParticleList
        from pytom.basic.files import read
        import os
        
        p = Particle('./testData/helicalSymmetryReference.em')
        pl = ParticleList()
        pl.append(p)
        
        vResult = read('./testData/helSymLeftResult.em')
        
        helSymObject = HelicalSymmetry(nfold = 4, symmetryAngle = 360/13.0, repeatShift = 21,isRightSymmetry = False)
        
        newPL = helSymObject.apply(pl)
        newPL.average('helSymTest.em',progressBar = False)
        
        tname = './helSymTest'
        v = read(tname+'.em')
        os.system('rm -f '+tname+'.em '+tname+'-PreWedge.em '+tname+'-WedgeSumUnscaled.em')
        
        assert vResult.equalsTo(v)
        

        
        
        
