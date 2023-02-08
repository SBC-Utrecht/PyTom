'''
Created on Jul 4, 2014

@author: foerster
'''

import unittest
import numpy as np
import pytom.voltools as vt
from pytom.basic.structures import Rotation
from pytom.angles.angleFnc import matToZXZ
from math import modf
from pytom.lib.pytom_numpy import vol2npy, npy2vol
from pytom.agnostic.tools import rotation_matrix, convert_angles, mat2ord
from scipy.ndimage import affine_transform
from pytom.agnostic.correlation import nxcc
from pytom.basic.transformations import rotate


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

        voltoolsmat1 = vt.utils.rotation_matrix(rotation=tuple(rot1.values()), rotation_order='rzxz')[:3, :3]
        voltoolsmat2 = vt.utils.rotation_matrix(rotation=tuple(rot2.values()), rotation_order='rzxz')[:3, :3]
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

        # This works for converting to relion
        zxz = np.array((9, 45, 114))
        zyz = convert_angles(zxz, rotation_order='zxz', return_order='zyz')
        zyz_relion = [a if pos else -a for a, pos in zip(zxz2zyz(*zxz), [0, 1, 0])]
        zyz_relion_alt = [-a for a in zyz]
        # zyz = mat2ord(rotation_matrix(zxz, rotation_order='zxz').T, return_order='zyz')
        # print(zyz_relion, zyz_relion_alt)
        relionmat = rotation_matrix(zyz_relion, rotation_order='zyz')
        relionmat_alt = rotation_matrix(zyz_relion_alt, rotation_order='zyz')
        # print(relionmat)
        # print(relionmat_alt)
        self.assertTrue(abs(np.array(zyz_relion) - np.array(zyz_relion_alt)).sum() < self.eps,
                        "something went wrong with forward backward transform to relion")
        self.assertTrue(np.abs(relionmat - relionmat_alt).sum() < self.eps,
                        "something went wrong with forward backward transform to relion")

        # I know this produces correct results in plotting neighbor density multiplying via simulation Vector
        # Vector uses dot(R, mat), the standard convention (I assume relion also does that)...
        vt_zyz_relion_custom = np.linalg.inv(vt.utils.rotation_matrix(rotation=[-a for a in zyz_relion],
                                                             rotation_order='rzyz')[:3, :3])
        # the rotation matrix is the same as the original pytom rotation matrix
        pytom_zxz = Rotation(z1=zxz[0], z2=zxz[2], x=zxz[1])
        pytom_zxz_mat = pytom_zxz.toMatrix()._matrix
        pytom_zxz_np = vol2npy(pytom_zxz_mat).copy(order='F')
        self.assertTrue(abs(vt_zyz_relion_custom - pytom_zxz_np).sum() < self.eps,
                        "something wrong with voltools conversion from relion versus pytom conversion")

        # but not the same as the voltools zxz matrix
        vt_zxz = vt.utils.rotation_matrix(rotation=zxz, rotation_order='rzxz')[:3, :3]
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
        # These are not equal
        # print(zyzToMat(zyz[0], zyz[2], zyz[1]))
        # print(pytom_zxz.toMatrix())
        # and this is not true ...
        # self.assertTrue(matToZXZ(zyzToMat(zyz[0], zyz[2], zyz[1])) == pytom_zxz,
        #                 "something went wrong with pytom internal angular conversion")

        vt_zyz = vt.utils.rotation_matrix(rotation=[-a for a in zyz], rotation_order='rzyz')[:3, :3]
        # if we take the negative values of the we get the correct matrix wit voltools
        # this means the rotation order is inverted z1, y, z2 == - (z2, y, z1)
        # print(np.linalg.inv(vt_zxz))
        # print(vt_zyz.T)
        self.assertTrue(np.abs(vt_zxz - vt_zyz).sum() < self.eps, "something wrong with rotation assumptions")

        # Angles dont have to be taken negative because pytom multiplies as dot(mat, R), while relion multiplies as
        # dot(R, mat). Relion also defines angles as from reference to particle, while pytom from particle to reference.
        # This is a tranpose operation but is taken care of by the order of multiplying z * y * z.

    def assertAnglesAlmostEqual(self, angles1, angles2, order):
        mat1 = rotation_matrix(angles1, order)
        mat2 = rotation_matrix(angles2, order)
        self.assertTrue(np.abs(mat1 - mat2).sum() < self.eps, "issue with angular conversion")

        # for i in [0, 1, 2]:
        #     self.assertAlmostEqual(angles1[i], angles2[i], places=3, msg="issue with angular conversions")

    def conversion(self, angles, order1, order2):
        # print('test post post')
        # test both pre and post multiplication
        new = convert_angles(angles, rotation_order=order1, return_order=order2, multiplication='post')
        res = convert_angles(new, rotation_order=order2, return_order=order1, multiplication='post')
        self.assertAnglesAlmostEqual(angles, res, order1)

        # print('test pre pre')
        new = convert_angles(angles, rotation_order=order1, return_order=order2, multiplication='pre')
        res = convert_angles(new, rotation_order=order2, return_order=order1, multiplication='pre')
        self.assertAnglesAlmostEqual(angles, res, order1)

        # print('test pre post')
        new = convert_angles(angles, rotation_order=order1, return_order=order2, multiplication='pre')
        res = convert_angles(new, rotation_order=order2, return_order=order1, multiplication='post')
        self.assertAnglesAlmostEqual(angles, res, order1)

        # print('test post pre')
        new = convert_angles(angles, rotation_order=order1, return_order=order2, multiplication='post')
        res = convert_angles(new, rotation_order=order2, return_order=order1, multiplication='pre')
        self.assertAnglesAlmostEqual(angles, res, order1)

    def object_rotation_test(self):
        # do numpy rotation
        obj = np.zeros((10, 11, 12))
        obj[4:5, 5:6, 7:8] += 10

        zxz = (150, 56, 210)

        zyz = convert_angles(zxz, rotation_order='zxz', return_order='zyz')

        zxz_post = rotation_matrix(zxz, rotation_order='zxz', multiplication='post')
        zyz_post = rotation_matrix(zyz, rotation_order='zyz', multiplication='post')

        obj_zxz = affine_transform(obj, zxz_post)
        obj_zyz = affine_transform(obj, zyz_post)

        self.assertTrue(np.abs(obj_zxz - obj_zyz).sum() < self.eps, "matrix from zxz or zyz are different")

        obj_back = affine_transform(obj_zxz, zxz_post.T)
        self.assertTrue(nxcc(obj, obj_back) > 0.8,
                        "forward and backward transform not the same")

        # compare with pytom basic
        obj_copy = obj.copy(order='F').astype(np.float32)
        # obj_pytom = vol(*obj.shape)
        obj_pytom = npy2vol(obj_copy, 3)
        obj_zxz_pytom = rotate(obj_pytom, zxz[0], x=zxz[1], z2=zxz[2])

        # back zxz angles
        zxz_back = mat2ord(zxz_post, return_order='zxz', multiplication='pre')
        obj_pytom_back = rotate(obj_zxz_pytom, zxz_back[0], x=zxz_back[1], z2=zxz_back[2])

        tmp = vol2npy(obj_pytom_back).copy(order='F')
        self.assertTrue(nxcc(tmp, obj) > 0.6, "forward backward zxz angles for pytom wrong")

    def angular_conversion_test(self):
        rot = (33, 20, 10)

        # test all forward and backward conversions
        options = ['xyz', 'xzy', 'yxz', 'yzx', 'zxy', 'zyx', 'xyx', 'xzx', 'yxy', 'yzy', 'zxz', 'zyz']
        conversions = []
        for oi in options:
            for oj in options:
                conversions.append((oi, oj))

        for c in conversions:
            self.conversion(rot, c[0], c[1])

    def voltools_agnostic(self):
        rot = (33, 20, 10)
        vtmat = vt.utils.rotation_matrix(rot, rotation_order='rzxz')[:3, :3]
        agmat = rotation_matrix(rot, rotation_order='zxz')
        self.assertTrue(np.abs(vtmat - agmat).sum() < self.eps, "incorrect assumption between pytom and voltools "
                                                                "angles")
        vtmat = vt.utils.rotation_matrix(rot, rotation_order='rzyz')[:3, :3]
        agmat = rotation_matrix(rot, rotation_order='zyz')
        # print(vtmat, agmat)
        self.assertTrue(np.abs(vtmat - agmat).sum() < self.eps, "incorrect assumption between pytom and voltools "
                                                                "angles")

    def angleFncTest(self):
        from pytom.angles.angleFnc import zyzToMat, matToZYZ, zxzToMat, matToZXZ
        order = [0, 2, 1]
        zyz = [9, 45, 114]
        zyz_conv = matToZYZ(zyzToMat(*[zyz[i] for i in order]))
        tmp1 = zyzToMat(*[zyz[i] for i in order])._matrix
        tmp2 = zyzToMat(*[zyz_conv[i] for i in order])._matrix
        zyzmat = vol2npy(tmp1).copy(order='F')
        zyzmat_conv = vol2npy(tmp2).copy(order='F')
        self.assertTrue(np.abs(zyzmat - zyzmat_conv).sum() < self.eps, "error in angles angleFnc angle conversion")

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
        self.angular_conversion_test()
        self.object_rotation_test()
        self.voltools_agnostic()
        # self.angleFncTest()


if __name__ == '__main__':
    unittest.main()


