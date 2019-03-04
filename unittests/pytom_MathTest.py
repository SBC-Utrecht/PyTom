#!/usr/bin/env pytom
'''
Created on Mar 26, 2010

@author: hrabe
'''

import unittest;

class pytom_MathTest(unittest.TestCase):
    def Matrix_Test(self):
        
        from pytom.tools.maths import Matrix;
        
        m = Matrix(3,3);
        
        assert m.getSizeX() == 3;
        assert m.getSizeY() == 3;
    
        #test __getitem__
        for i in range(9):
            assert m[int(i//3),int(i%3)] == 0;
            
        #test __setitem__
        for i in range(9):
            m[int(i//3),int(i%3)] = i+1;
            assert m[int(i//3),int(i%3)] == i+1;
            
        #test rows and columns
        assert [1,2,3] == m.getRow(0);
        assert [1,4,7] == m.getColumn(0);
        
        assert [4,5,6] == m.getRow(1);
        assert [2,5,8] == m.getColumn(1);
        
        assert [7,8,9] == m.getRow(2);
        assert [3,6,9] == m.getColumn(2);
        
                
        m2 = m * m;
        assert [30.0, 36.0, 42.0] == m2.getRow(0);
        assert [66.0, 81.0, 96.0] == m2.getRow(1);
        assert [102.0, 126.0, 150.0] == m2.getRow(2);
        
                
    def Identity_Test(self):
        
        from pytom.tools.maths import Identity;
        
        
        i = Identity(3,3);
        i2 = Identity(3,3);
        
        assert (i*i2) == i;
        assert i.isIdentity();
        
        
    def RotationMatrix_Test(self):
        from pytom.tools.maths import XRotationMatrix,YRotationMatrix,ZRotationMatrix,Identity;
        
        i = Identity(3,3);
        x = XRotationMatrix(0);
        y = YRotationMatrix(0);
        z = ZRotationMatrix(0);
        
        assert x == i and x.isIdentity();
        assert y == i and y.isIdentity();
        assert z == i and z.isIdentity();
    
    def Pcacov_Test(self):
        from pytom_volume import vol
        
        v = vol(4,4,1)
        v.setV(7,0,0,0)
        v.setV(1,1,0,0)
        v.setV(11,2,0,0)
        v.setV(11,3,0,0)
        
        v.setV(26,0,1,0)
        v.setV(29,1,1,0)
        v.setV(56,2,1,0)
        v.setV(31,3,1,0)
        
        v.setV(6,0,2,0)
        v.setV(15,1,2,0)
        v.setV(8,2,2,0)
        v.setV(8,3,2,0)
        
        v.setV(60,0,3,0)
        v.setV(52,1,3,0)
        v.setV(20,2,3,0)
        v.setV(47,3,3,0)
        
        from pytom.tools.maths import pcacov
        
        [coeff,latent,explained] = pcacov(v)
        
        from numpy import zeros,sum
        
        z = zeros(4)
        
        z[0] = 70.06892395
        z[1] = 22.51976204
        z[2] = 5.57571507
        z[3] = 1.83560765
        
        assert sum(explained == z) < 0.000000001
        
    def sumRotations_T(self,r1,r2,result):
        
        assert r1*r2 == result
        
    def sumRotations_Test(self):
        from pytom.basic.structures import Rotation
        self.sumRotations_T(Rotation(0,0,0),Rotation(0,0,0),Rotation(0.0,0.0,0.0))
        self.sumRotations_T(Rotation(0,0,0),Rotation(10,10,10),Rotation(10,10,10))
        self.sumRotations_T(Rotation(10,0,0),Rotation(0,0,10),Rotation(0,10,10))
        self.sumRotations_T(Rotation(10,20,30),Rotation(0,0,10),Rotation(7.78210641642,22.6953896502,39.8822808113))
        
        