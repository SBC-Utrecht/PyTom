import unittest

class pytom_FRMTest(unittest.TestCase):

    def FRM_Test(self):
        import swig_frm
        from sh_alignment.frm import frm_align
        from pytom_volume import vol, rotate, shift
        from pytom.basic.structures import Rotation, Shift

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

        print(Rotation(ang), rotation)

        assert Rotation(ang) == rotation
        assert Shift(pos) == shiftO
        
