'''
Created on Jul 4, 2014

@author: foerster
'''

import unittest
import numpy as np
from pytom.basic.structures import Rotation
from pytom.angles.angleFnc import matToZXZ
from math import modf
from pytom.voltools.utils import rotation_matrix
from pytom_numpy import vol2npy, npy2vol


class pytom_RotationTest(unittest.TestCase):

    def setUp(self):
        """set up"""
        #self.ang1=Rotation(z1=-016,z2=-114,x=0.2)
        self.Pole1 = Rotation(z1=9,z2=-114,x=0.)
        self.Pole2 = Rotation(z1=9,z2=-114,x=180.)
        self.ang1=Rotation(z1=-16,z2=-114,x=-180)
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

    def voltools_pytom_test(self):
        """
        Check the relation between pytom and voltools rotation matrices.
        We want to rotate something first with rot1 and then with rot2.
        If coordinates are rotated with:
        - dot(rot, R), then we need to do dot(rot2, rot1)  => conv1
        - dot(R, rot), then we need to do dot(rot1, rot2)  => conv2
        We can find the link through rot_conv1 = linalg.inv(rot_conv2)
        """
        rot1 = {'z1': 9,
                'x': 45,
                'z2': 114}
        rot2 = {'z1': 56,
                'x': 150,
                'z2': 30}

        pytomrot1 = Rotation(z1=rot1['z1'], z2=rot1['z2'], x=rot1['x'])
        pytomrot2 = Rotation(z1=rot2['z1'], z2=rot2['z2'], x=rot2['x'])
        pytom_final = pytomrot2 * pytomrot1
        pytom_final_matrix = pytom_final.toMatrix()._matrix
        pytom_final_matrix_np = vol2npy(pytom_final_matrix).copy(order='F')

        voltoolsmat1 = rotation_matrix(rotation=tuple(rot1.values()), rotation_order='rzxz')[:3, :3]
        voltoolsmat2 = rotation_matrix(rotation=tuple(rot2.values()), rotation_order='rzxz')[:3, :3]
        voltools_final = np.dot(voltoolsmat1, voltoolsmat2)
        voltools_final_compat = np.linalg.inv(voltools_final)

        # these are already inverted based on the order of zxz multiplication
        # print(pytomrot1.toMatrix())
        # print(voltoolsmat1)

        # print(pytom_final_matrix_np)
        # print(voltools_final_compat)
        diff_sum = np.abs(pytom_final_matrix_np - voltools_final_compat).sum()
        self.assertTrue(diff_sum < self.eps, "conversion between pytom and voltools rotation matrices failed")

    def conversion_to_relion(self):
        from pytom.agnostic.tools import zxz2zyz, zyz2zxz
        from pytom.angles.angleFnc import matToZYZ, zyzToMat, matToZXZ

        zxz = np.array((9, 45, 114))
        tmp = zxz2zyz(*zxz)
        zyz_relion = (- tmp[0], tmp[1], - tmp[2])
        self.assertTrue(abs(np.array(zyz2zxz(- zyz_relion[0], zyz_relion[1], - zyz_relion[2])) - zxz).sum() < self.eps,
                        "something went wrong with forward backward transform to relion")

        # I know this produces correct results in plotting neighbor density multiplying via simulation Vector
        # Vector uses dot(R, mat), the standard convention (I assume relion also does that)...
        vt_zyz_relion_custom = np.linalg.inv(rotation_matrix(rotation=[-a for a in zyz_relion],
                                                             rotation_order='rzyz')[:3, :3])
        # the rotation matrix is the same as the original pytom rotation matrix
        pytom_zxz = Rotation(z1=zxz[0], z2=zxz[2], x=zxz[1])
        pytom_zxz_mat = pytom_zxz.toMatrix()._matrix
        pytom_zxz_np = vol2npy(pytom_zxz_mat).copy(order='F')
        self.assertTrue(abs(vt_zyz_relion_custom - pytom_zxz_np).sum() < self.eps,
                        "something wrong with voltools conversion from relion versus pytom conversion")

        # but not the same as the voltools zxz matrix
        vt_zxz = rotation_matrix(rotation=zxz, rotation_order='rzxz')[:3, :3]
        # vt_zyz_relion = rotation_matrix(rotation=zyz_relion, rotation_order='rzyz')[:3, :3]
        # print(pytom_zxz_np)
        # print(vt_zxz)
        # print(np.linalg.inv(vt_zyz_relion_custom))

        # print(vt_zyz_relion)

        # This produces the same results as using the zxz2zyz from agnostic with the z angle inversion
        tmp = matToZYZ(pytom_zxz.toMatrix())
        zyz = (tmp[0], tmp[1], tmp[2])
        self.assertTrue(np.abs(np.array(zyz) - np.array(zyz_relion)).sum() < self.eps,
                        "something wrong with conversion assumptions")

        vt_zyz = rotation_matrix(rotation=[-a for a in zyz], rotation_order='rzyz')[:3, :3]
        # if we take the negative values of the we get the correct matrix wit voltools
        # this means the rotation order is inverted z1, y, z2 == - (z2, y, z1)
        # print(np.linalg.inv(vt_zxz))
        # print(vt_zyz.T)
        self.assertTrue(np.abs(vt_zxz - vt_zyz).sum() < self.eps, "something wrong with rotation assumptions")

        # Angles dont have to be taken negative because pytom multiplies as dot(mat, R), while relion multiplies as
        # dot(R, mat). Relion also defines angles as from reference to particle, while pytom from particle to reference.
        # This is a tranpose operation but is taken care of by the order of multiplying z * y * z.

    def runTest(self):
        """
        """
        self.multiplicationTest()
        self.matTestQ1()
        self.matTestQ2()
        self.matTestQ3()
        self.matTestQ4()
        self.matPoleTest1()
        self.matPoleTest2()
        self.voltools_pytom_test()
        self.conversion_to_relion()


if __name__ == '__main__':
    unittest.main()


