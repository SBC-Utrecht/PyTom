#!/usr/bin/env pytom
'''
Created on Mar 26, 2010
redone Jan 8, FF

@author: hrabe, FF
'''

import unittest;

class pytom_MathTest(unittest.TestCase):
    def test_Matrix(self):
        
        from pytom.tools.maths import Matrix;
        
        m = Matrix(3,3);
        
        self.assertTrue( m.getSizeX() == 3,
            msg='Matrix dim incorrect')
        self.assertTrue( m.getSizeY() == 3,
            msg='Matrix dim incorrect')
    
        #test __getitem__
        for i in range(9):
            assert m[int(i/3),int(i%3)] == 0;
            
        #test __setitem__
        for i in range(9):
            m[int(i/3),int(i%3)] = i+1;
            assert m[int(i/3),int(i%3)] == i+1;
            
        #test rows and columns
        self.assertTrue( [1,2,3] == m.getRow(0),
            msg='rows incorrect')
        self.assertTrue( [1,4,7] == m.getColumn(0),
            msg='rows incorrect')
        self.assertTrue( [4,5,6] == m.getRow(1),
            msg='rows incorrect')
        self.assertTrue( [2,5,8] == m.getColumn(1),
            msg='rows incorrect')
        self.assertTrue( [7,8,9] == m.getRow(2),
            msg='rows incorrect')
        self.assertTrue( [3,6,9] == m.getColumn(2),
            msg='rows incorrect')
                
        m2 = m * m;
        self.assertTrue( [30.0, 36.0, 42.0] == m2.getRow(0),
            msg='multi incorrect')
        self.assertTrue( [66.0, 81.0, 96.0] == m2.getRow(1),
            msg='multi incorrect')
        self.assertTrue( [102.0, 126.0, 150.0] == m2.getRow(2),
            msg='multi incorrect')
        
                
    def test_Identity(self):
        from pytom.tools.maths import Identity;
        
        i = Identity(3,3);
        i2 = Identity(3,3);
        
        self.assertTrue( (i*i2) == i, msg='Identity incorrect')
        self.assertTrue( i.isIdentity(), msg='Identity incorrect')
        
        
    def test_RotationMatrix(self):
        from pytom.tools.maths import XRotationMatrix,YRotationMatrix,ZRotationMatrix,Identity;
        
        i = Identity(3,3);
        x = XRotationMatrix(0);
        y = YRotationMatrix(0);
        z = ZRotationMatrix(0);
        
        self.assertTrue( x == i and x.isIdentity(),
                    msg='Rotation incorrect')
        self.assertTrue( y == i and y.isIdentity(),
                    msg='Rotation incorrect')
        self.assertTrue( z == i and z.isIdentity(),
                    msg='Rotation incorrect')
    
    def test_Pcacov(self):
        from pytom.lib.pytom_volume import vol
        from pytom.tools.maths import pcacov
        from numpy import zeros,sum
        
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
        
        [coeff,latent,explained] = pcacov(v)
        
        z = zeros(4)
        
        z[0] = 70.06892395
        z[1] = 22.51976204
        z[2] = 5.57571507
        z[3] = 1.83560765
        
        self.assertTrue( sum(explained == z) < 0.000000001,
                    msg='PCA-Cov incorrect')
        
    def sumRotations_T(self,r1,r2,result):
        self.assertTrue( r1*r2 == result, msg='Rotation incorrect')
        
    def test_sumRotations(self):
        from pytom.basic.structures import Rotation

        self.sumRotations_T(Rotation(0,0,0),Rotation(0,0,0),Rotation(0.0,0.0,0.0))
        self.sumRotations_T(Rotation(0,0,0),Rotation(10,10,10),Rotation(10,10,10))
        self.sumRotations_T(Rotation(10,0,0),Rotation(0,0,10),Rotation(0,10,10))
        self.sumRotations_T(Rotation(10,20,30),Rotation(0,0,10),Rotation(7.78210641642,22.6953896502,39.8822808113))
        
        
if __name__ == '__main__':
    unittest.main()
