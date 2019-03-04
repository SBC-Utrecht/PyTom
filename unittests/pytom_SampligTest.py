'''
Created on Jul 7, 2014

@author: foerster
'''

import unittest
import os
from pytom.alignment.alignmentFunctions import bestAlignment
from pytom.basic.structures import Wedge, Particle
from pytom.angles.localSampling import LocalSampling
from pytom.basic.files import read
from pytom.tools.maths import rotation_distance
from pytom.score.score import nxcfScore
from pytom_volume import vol, transform# as transform #developers: your can also import transformSpline for more accurate rotation!

class pytom_SamplingTest(unittest.TestCase):

    def setUp(self):
        """set up"""
	self.testFile = 'testData/emd_1480.map.em_bin_4.em'
	self.ref = read(self.testFile)
	self.particle = read(self.testFile)
	self.rotAngles = []
	self.rotAngles.append([0.,0.,0.])
	self.rotAngles.append([10.,20.,30.])
	self.rotAngles.append([0.,0.,90.])
	self.rotAngles.append([90.,270.,90.])
	self.startAngles = []
	self.startAngles.append([10.,20.,30.])
	self.startAngles.append([30.,40.,10.])
	self.startAngles.append([350.,20.,110.])
	self.startAngles.append([70.,260.,75.])
        self.scoreObj = nxcfScore(-1000)
	self.eps = 0.001
	self.bestRot = [10,20,30]
	self.iniRot = [0,0,0]
	self.ccc = -1000.
	self.nshells = 3
	self.increments = [7.5, 5., 2., 1., .5, .25, .15, .1]

    def match(self, nshells=1, increment=1., z1Start=0., z2Start=0., xStart=0.):
        """
	
        """
        rotations = LocalSampling(shells=self.nshells,increment=increment,
	     z1Start=z1Start, z2Start=z2Start, xStart=xStart)
	#rotations.printRotations()
        bestPeak =bestAlignment(particle=self.particle, reference=self.ref, 
	    referenceWeighting="None",wedgeInfo=Wedge(),rotations=rotations,
            scoreObject=self.scoreObj, mask=None,preprocessing=None,progressBar=False,binning=1,
	    bestPeak = None,verbose = False)
	return bestPeak

    def sampleLocally(self, rotAngle, startAngle):
        """
	iterate with decreasing angular increments
	"""
	self.particle = vol( self.ref.sizeX(), self.ref.sizeY(), self.ref.sizeZ())
	transform( self.ref, self.particle, rotAngle[0], rotAngle[1], rotAngle[2], 
	    int(self.ref.sizeX()/2), int(self.ref.sizeX()/2), 
	    int(self.ref.sizeX()/2), 0,0,0,0,0,0)
	self.particle.write('xxx.em')
	self.bestRot = [startAngle[0], startAngle[1], startAngle[2]]
	for (ii, increment) in enumerate(self.increments):
	    bestPeak = self.match(nshells=self.nshells, increment=increment, 
	        z1Start=self.bestRot[0], z2Start=self.bestRot[1], xStart=self.bestRot[2])
	    self.bestRot = bestPeak.getRotation()
	    self.ccc = bestPeak.getScoreValue()
	    self.trans = bestPeak.getShift()
	    #print bestPeak
	return bestPeak



    def eval(self):
        """
	"""
	for (ii, rotAngle) in enumerate(self.rotAngles):
	    #print "\n\n--------------------------------------"
	    #print rotAngle
	    #print "--------------------------------------"
	    startAngle = self.startAngles[ii]
	    self.sampleLocally(rotAngle, startAngle)
	    dist = rotation_distance(self.bestRot, rotAngle)
	    #print "rotAngle ="+str(rotAngle)+": "+str(dist)
            self.assertTrue( abs(dist) < 0.2, 'Accuracy of Local Sampling > 0.2 deg')

        
    def runTest(self):
        self.eval()

if __name__ == '__main__':
    unittest.main()


