import unittest


class pytom_2DImageTest(unittest.TestCase):
    def setUp(self):

        from pytom.basic.functions import initSphere
        from pytom.reconstruction.imageStructures import ImageStack, Image
        from pytom.basic.transformations import general_transform2d
        from pytom.simulation.support import add_white_noise
        import random

        random.seed(0)
        self.dim = 32
        self.nim = 41
        self.snr= .1
        self.radius = 7.

        self.sph = initSphere(size_x=self.dim, size_y=self.dim, size_z=1, 
            radius=self.radius, smooth=1., 
            maxradius=0, cent=None, filename='')
        self.imageStack = ImageStack()

        # random shift
        self.shiftX=[]
        self.shiftY=[]
        sumX=0.
        sumY=0.
        for ii in range(0,self.nim):
            self.shiftX.append(random.gauss(mu=0., sigma=2.))
            self.shiftY.append(random.gauss(mu=0., sigma=2.))
        sumX = sumX/self.nim
        sumY = sumY/self.nim
        for ii in range(0,self.nim):
            self.shiftX[ii] = self.shiftX[ii]-sumX
            self.shiftY[ii] = self.shiftY[ii]-sumY
        
        # apply shifts
        for ii in range(0,self.nim):
            image = Image( filename=None, boxCoords=[0,0,0], dims=[self.dim,self.dim,1],
                    shiftX=0, shiftY=0, appliedShiftX=0, appliedShiftY=0, rotation=0,
            index = ii, verbose=False)
            tmp = general_transform2d( v=self.sph, rot=0., 
                    shift=[self.shiftX[ii], self.shiftY[ii]], 
                    scale=1., order=[2, 1, 0], crop=True)
            image.data = add_white_noise(volume=tmp,SNR=self.snr)
            self.imageStack.append(image=image)
    
        self.imageStack.normalize( normtype="StdMean")
        self.imageStack.bandpass( lowfreq=1., hifreq=6., smooth=2., bpf=None)
        # smoothen edges
        self.imageStack.taper_edges( width=self.dim/8.)

        self.imageStack.exMaxAlign( niter=10, mask=None)
        self.imageStack.write('input.em')
        self.imageStack.writeWorkingCopies('output.em')


    def eval(self):
        """
        compare determined shifts to ground truth
        """
        import os
        from math import sqrt
        rms = 0.
        for ii in range(0,self.nim):
            tline = "inputXY=%5.2f" %self.shiftX[ii]
            tline = tline+", %5.2f" %self.shiftY[ii]
            tline = tline+"; outputXY=%5.2f" %self.imageStack.images[ii].shiftX
            tline = tline+", %5.2f" %self.imageStack.images[ii].shiftY
            #print tline
            tmp = (self.shiftX[ii]-self.imageStack.images[ii].shiftX)**2
            tmp = (self.shiftY[ii]-self.imageStack.images[ii].shiftY)**2 + tmp
            rms = rms + tmp
        rms = sqrt(rms / self.nim)
        tline = "RMSD = %5.2f" %rms
        print(tline)
        self.assertTrue( rms < 1.0, 'RMSD of ExMax alignment too large!')

        os.system('rm input.em output.em')

    def runTest(self):
        self.eval()

if __name__ == '__main__':
    unittest.main()



