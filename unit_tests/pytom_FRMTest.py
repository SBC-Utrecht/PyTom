import unittest

class pytom_FRMTest(unittest.TestCase):

    def mySetUp(self):
        from helper_functions import create_RandomParticleList
        print("running setup")
        self.reffile = './testData/ribo.em'
        self.pl_filename = 'pl.xml'
        self.pdir = './testparticles' 
        self.pl = create_RandomParticleList( reffile=self.reffile, pl_filename=self.pl_filename, 
                  pdir=self.pdir, nparticles=10)

        # set parameters
        self.settings = {}
        self.settings["frequency"] = 10
        self.settings["peak_offset"] = 5
        self.settings["binning"] = 2
        self.settings["niteration"] = 2
        self.settings["mask"] = './testData/ribo_mask.em'
        self.settings["pixsize"] = 2.7
        self.settings["adres"]   = 0.00
        self.settings["rescrit"] = 0.5
        self.settings["jobfile"] = 'xxx.xml'

    #def cleanUp(self):
    def myCleanUp(self):
        """
        check that files are written and remove them
        """
        from helper_functions import cleanUp_RandomParticleList
        from os import system

        print("running cleanup")
        for ii in range(0,self.settings["niteration"]):
            tline=('rm average_iter'+str(ii)+'*.em')
            system(tline)
            tline=('rm fsc_'+str(ii)+'_odd.em')
            system(tline)
            tline=('rm fsc_'+str(ii)+'_even.em')
            system(tline)
            tline=('rm aligned_pl_iter'+str(ii)+'.xml')
            system(tline)
        system('rm '+self.settings["jobfile"])
        #cleanUp_RandomParticleList( pl_filename=self.pl_filename, pdir=self.pdir)
        #self.remove_file( filename=self.settings["jobfile"])

    def remove_file(self, filename):
        """
        assert that file exists end remove it
        """
        from os import remove
        from os import path

        filecheck = path.exists(filename)
        self.assertTrue( filecheck, msg="file "+filename+" does not exist")
        if filecheck:
            remove(filename)

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

    def test_commandLine(self):
        """test that script works"""
        from os import system
        self.mySetUp()
        cmd = (r"printf '0\n"+
                     self.pl_filename+r'\n'+
                     self.reffile+r'\n'+
                     self.settings["mask"]+r'\n'+
                     str(self.settings["frequency"])+r'\n'+
                     str(self.settings["peak_offset"])+r'\n'+
                     str(self.settings["niteration"])+r'\n'+
                     str(self.settings["binning"])+r'\n'+
                     str(self.settings["pixsize"])+r'\n'+
                     str(self.settings["adres"])+r'\n'+
                     str(self.settings["rescrit"])+r'\n'+
                     self.settings["jobfile"] + "' | ../bin/pytom ../frm/createJob.py")
        e = system(cmd)
        cmd = 'mpirun -np 2 ../bin/pytom ../frm/FRMAlignment.py -j '+self.settings["jobfile"]
        e = system(cmd)
        myCleanUp()
        #try:
        #    e = system(cmd)
        #except:
        #    print('hello')
        #self.assertTrue( e == 256, 'no error raised!')

        
if __name__ == '__main__':
    unittest.main()

