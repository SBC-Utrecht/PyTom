import unittest
from pytom_volume import vol
from numpy import array
from pytom.basic.transformations import general_transform2d
from pytom.tools.maths import rotate_vector


class pytom_TransformationTest(unittest.TestCase):
    def setUp(self):
        """set up"""
        # a point somewhere
        self.myVol = vol(256,256,1)
        self.myVol.setAll(0.)
        self.myVol.setV(1.,200,128,0)
        self.r = array([200-128,0])
        self.origin = vol(256,256,1)
        self.origin.setAll(0.)
        self.origin.setV(1.,128,128,0)
        self.m = 0.8

    def test_Origin(self):
        """
        test that central voxel remains invariant upon magification
        """
        newvol = general_transform2d(self.origin, rot=0., shift=None, scale=self.m)
        val = newvol.getV(128,128,0)
        self.assertTrue( val> 0.01, 'origin changes upon magnification change!')

    def test_Transform(self):
        """
        """
        # first magnify by m, then rotate by alpha
        alpha = 25.
        self.r_trans = rotate_vector(self.r, alpha)
        self.r_trans = self.r_trans * self.m
        r_trans_round = array([0.,0.])
        tline = "Approximate position after rotation and magnification: "
        for ii in range(0,2):
            r_trans_round[ii] = int(round(self.r_trans[ii]) + 128)
            tline = tline + ("%6d " %r_trans_round[ii])
        print(tline)
        newvol = general_transform2d(self.myVol, rot=alpha, shift=None, scale=self.m)
        val = newvol.getV(int(r_trans_round[0]),int(r_trans_round[1]),0)
        print("value at approximate interpolated position: %5.3f" %val)
        tline = ("next neighbors: %5.3f, " %newvol.getV(int(r_trans_round[0])-1,int(r_trans_round[1]),0))
        tline = tline + ("%5.3f, " %newvol.getV(int(r_trans_round[0])-1,int(r_trans_round[1]),0))
        tline = tline + ("%5.3f, " %newvol.getV(int(r_trans_round[0]),int(r_trans_round[1])-1,0))
        tline = tline + ("%5.3f, " %newvol.getV(int(r_trans_round[0])+1,int(r_trans_round[1]),0))
        tline = tline + ("%5.3f, " %newvol.getV(int(r_trans_round[0]),int(r_trans_round[1])+1,0))
        tline = tline + ("%5.3f, " %newvol.getV(int(r_trans_round[0]+1),int(r_trans_round[1])+1,0))
        tline = tline + ("%5.3f, " %newvol.getV(int(r_trans_round[0]+1),int(r_trans_round[1])-1,0))
        tline = tline + ("%5.3f, " %newvol.getV(int(r_trans_round[0]-1),int(r_trans_round[1])+1,0))
        tline = tline + ("%5.3f\n" %newvol.getV(int(r_trans_round[0]-1),int(r_trans_round[1])-1,0))
        print(tline)
        self.assertTrue( val> 0.01, 'wrong point !=0!')

    def test_resize2D(self):
        """
        test re-sizing in Fourier space
        """
        from pytom.basic.transformations import resize
        from pytom_volume import vol
        from pytom.basic.fourier import fft

        dim = 32
        px = 11
        py = 19
        scf = dim*dim
        myVol = vol(dim,dim,1)
        myVol.setAll(0.)
        myVol.setV(1., px, py, 0)
        #fmyVol = fft(myVol)
        (resizeVol, resizefVol) = resize(volume=myVol, factor=2., interpolation='Fourier')
        #resizeVol.write('test1.em')
        ftresizeVol = fft(data=resizeVol)
        for ix in range(resizefVol.sizeX()):
            for iy in range(resizefVol.sizeY()):
                diff = ftresizeVol.getV(ix,iy,0) - scf*4*resizefVol.getV(ix,iy,0)
                self.assertTrue(expr=abs(diff) < .05, msg="inconsistency FFT/IFFT for magnification")
        (resizeVol, resizefVol) = resize(volume=resizeVol, factor=.5, interpolation='Fourier')
        from pytom_volume import variance
        diff = myVol - resizeVol
        self.assertTrue(expr=variance(diff, False) < .0000001, msg="2D image before and after rescales differs")


    def test_resize3D(self):
        """
        test 3D re-sizing in Fourier space
        """
        from pytom.basic.transformations import resize
        from pytom_volume import vol

        dim = 32
        px = 11
        py = 19
        pz = 21
        myVol = vol(dim,dim,dim)
        myVol.setAll(0.)
        myVol.setV(1., px, py, pz)
        (resizeVol, resizefVol) = resize(volume=myVol, factor=2., interpolation='Fourier')
        (resizeVol, resizefVol) = resize(volume=resizeVol, factor=.5, interpolation='Fourier')
        from pytom_volume import variance
        diff = myVol - resizeVol
        self.assertTrue(expr=variance(diff, False) < .0000001, msg="2D image before and after rescales differs")

    def runTest(self):
        self.test_Origin()
        self.test_Transform()
        self.test_resize2D()
        self.test_resize3D()
        self.test_FullReducedTrafo()


if __name__ == '__main__':
    unittest.main()

