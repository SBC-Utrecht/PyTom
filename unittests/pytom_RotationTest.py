'''
Created on Jul 4, 2014

@author: foerster
'''

import unittest
from pytom.basic.structures import Rotation
from pytom.angles.angleFnc import matToZXZ
from math import modf


class pytom_RotationTest(unittest.TestCase):

    def setUp(self):
        """set up"""
	#self.ang1=Rotation(z1=-016,z2=-114,x=0.2)
	self.Pole1 = Rotation(z1=9,z2=-114,x=0.)
	self.Pole2 = Rotation(z1=9,z2=-114,x=180.)
	self.ang1=Rotation(z1=-016,z2=-114,x=-180)
	self.ang2=Rotation(z1=0,z2=5.,x=0.0)
	self.eps = .001

    def matTestQ1(self):
        """
	test that matToZXZ works for 0<z2<180
        """
	z1=9
	z2=170.
	x=10.
	self.rot1 = Rotation(z1=z1,z2=z2,x=x)
	self.rot2 = Rotation(z1=0.,z2=0.,x=0.)
	newrot    = self.rot2*self.rot1
        [nz1, nz2, nx] = matToZXZ(newrot.toMatrix(), inRad=False) 
        self.assertTrue( abs(nx -x) < self.eps, 'matTestQ1: nx and x differ')
        self.assertTrue( abs(nz1 -z1) < self.eps, 'matTestQ1: nz1 and z1 differ')
        self.assertTrue( abs(nz2 -z2) < self.eps, 'matTestQ1: nz2 and z2 differ')

    def matTestQ2(self):
        """
	test that matToZXZ works for 180<z1<360
        """
	z1=189.
	z2=170.
	x=10.
	self.rot1 = Rotation(z1=z1,z2=z2,x=x)
	self.rot2 = Rotation(z1=0.,z2=0.,x=0.)
	newrot    = self.rot2*self.rot1
        [nz1, nz2, nx] = matToZXZ(newrot.toMatrix(), inRad=False) 
        self.assertTrue( abs(nx -x) < self.eps, 'matTestQ2: nx and x differ')
        self.assertTrue( abs(nz1 -z1) < self.eps, 'matTestQ2: nz1 and z1 differ')
        self.assertTrue( abs(nz2 -z2) < self.eps, 'matTestQ2: nz2 and z2 differ')

    def matTestQ3(self):
        """
	test that matToZXZ works for 180<z2<360
        """
	z1=189.
	z2=190.
	x=10.
	self.rot1 = Rotation(z1=z1,z2=z2,x=x)
	self.rot2 = Rotation(z1=0.,z2=0.,x=0.)
	newrot    = self.rot2*self.rot1
        [nz1, nz2, nx] = matToZXZ(newrot.toMatrix(), inRad=False) 
        self.assertTrue( abs(nx -x) < self.eps, 'matTestQ3: nx and x differ')
        self.assertTrue( abs(nz1 -z1) < self.eps, 'matTestQ3: nz1 and z1 differ')
        self.assertTrue( abs(nz2 -z2) < self.eps, 'matTestQ3: nz2 and z2 differ')

    def matTestQ4(self):
        """
	test that matToZXZ works for 180<z2<360
        """
	z1=169.
	z2=190.
	x=10.
	self.rot1 = Rotation(z1=z1,z2=z2,x=x)
	self.rot2 = Rotation(z1=0.,z2=0.,x=0.)
	newrot    = self.rot2*self.rot1
        [nz1, nz2, nx] = matToZXZ(newrot.toMatrix(), inRad=False) 
        self.assertTrue( abs(nx -x) < self.eps, 'matTestQ3: nx and x differ')
        self.assertTrue( abs(nz1 -z1) < self.eps, 'matTestQ3: nz1 and z1 differ')
        self.assertTrue( abs(nz2 -z2) < self.eps, 'matTestQ3: nz2 and z2 differ')


    def matPoleTest1(self):
        """
	check that conversion to angle from matrix works for x == 0
	"""
	z1=9
	z2=360.-114
	x=0.
	self.rot1 = Rotation(z1=z1,z2=z2,x=x)
	self.rot2 = Rotation(z1=0.,z2=0.,x=0.)
	newrot    = self.rot2*self.rot1
        [nz1, nz2, nx] = matToZXZ(newrot.toMatrix(), inRad=False) 
        self.assertTrue( nx == 0.0, 'Pole Test 1 failed: wrong x-angle')
        self.assertTrue( abs(nz2 + nz1 - z2 -z1) < 0.00001, 
	    'Pole Test 2 failed: wrong assignment z1 and z2')

    def matPoleTest2(self):
        """
	check that conversion to angle from matrix works for x == 180
	"""
	z1=9
	z2=-114
	x=180.
	self.rot1 = Rotation(z1=z1,z2=z2,x=x)
	self.rot2 = Rotation(z1=0.,z2=0.,x=0.)
	newrot    = self.rot2*self.rot1
        [nz1, nz2, nx] = matToZXZ(newrot.toMatrix(), inRad=False) 
        self.assertTrue( nx == 180.0, 'Pole Test 2 failed: wrong x-angle')
        self.assertTrue( abs(modf((nz2 - nz1 + 360.)/360.)[0]*360 - 
	    modf((z2 -z1 + 360.)/360.)[0]*360) < 0.00001, 
	    'Pole Test 2 failed: z2 z1 not correct')

        
    def multiplicationTest(self):
        """
        """
	multa1a2=self.ang2*self.ang1


    def eval(self):
        """
	"""
        self.multiplicationTest()
        self.matTestQ1()
        self.matTestQ2()
        self.matTestQ3()
        self.matTestQ4()
	self.matPoleTest1()
	self.matPoleTest2()

        
    def runTest(self):
        self.eval()

if __name__ == '__main__':
    unittest.main()


