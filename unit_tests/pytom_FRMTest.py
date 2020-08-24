import unittest

class pytom_FRMTest(unittest.TestCase):

    def test_FRM(self):
        import swig_frm
        from sh_alignment.frm import frm_align
        from pytom_volume import vol, rotate, shift
        from pytom.basic.structures import Rotation, Shift
        from pytom.tools.maths import rotation_distance

        v = vol(32,32,32)
        v.setAll(0)
        vMod = vol(32,32,32)
        vRot = vol(32,32,32)
        v.setV(1,10,10,10)
        v.setV(1,20,20,20)
        v.setV(1,15,15,15)
        v.setV(1,7,21,7)
        
        rotation = Rotation(10,20,30)
        shiftO = Shift(1,-3,5)
        
        rotate(v,vRot,rotation.getPhi(),rotation.getPsi(),rotation.getTheta())
        shift(vRot,vMod,shiftO.getX(),shiftO.getY(),shiftO.getZ())
        
        pos, ang, score = frm_align(vMod, None, v, None, [4, 64], 10)
        rotdist = rotation_distance(ang1=rotation, ang2=ang)
        diffx = shiftO[0] - (pos[0] - 16)
        diffy = shiftO[1] - (pos[1] - 16)
        diffz = shiftO[2] - (pos[2] - 16)

        self.assertTrue( rotdist < 5., msg='Rotations are different')
        self.assertTrue( diffx < .5, msg='x-difference > .5')
        self.assertTrue( diffy < .5, msg='y-difference > .5')
        self.assertTrue( diffz < .5, msg='z-difference > .5')
        
if __name__ == '__main__':
    unittest.main()
