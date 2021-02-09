"""
test a number of filters, also for reconstruction
"""
import unittest

class pytom_FilterTest(unittest.TestCase):
    
    def test_wedgeFilter(self):
        """
        test wedge filter functions and their consistency
        """
        from pytom_volume import vol
        from pytom.basic.structures import WedgeInfo,Rotation,Shift
        from pytom.basic.filter import filter
        from pytom.basic.fourier import fft, ifft
        from pytom_volume import complexRealMult

        v = vol(32,32,32)
        v.setAll(0.0)
        v.setV(1.,16,16,16)

        wedgeInfo = WedgeInfo(30.0)#,[10.0,20.0,30.0],0.0)
        wvol = wedgeInfo.returnWedgeVolume(v.sizeX(), v.sizeY(), v.sizeZ())
        wfil = wedgeInfo.returnWedgeFilter(v.sizeX(), v.sizeY(), v.sizeZ())
        vfil1 = wedgeInfo.apply(v)
        vfil2 = filter(v,wfil)

        fv = fft(v)
        fvol3 = complexRealMult(fv,wvol)
        vfil3 = ifft(fvol3)
        vfil3.shiftscale(0.0,1/float(v.sizeX()*v.sizeY()*v.sizeZ()))

        self.assertAlmostEqual(vfil1.getV(16,16,16),vfil2[0].getV(16,16,16),2,"Wedge Filter Inconsistency")
        self.assertAlmostEqual(vfil1.getV(16,16,16),vfil3.getV(16,16,16),2,"Wedge Filter Inconsistency")

    def test_wedgeRotation(self):
        """
        """
        from pytom_volume import vol
        from pytom.basic.structures import Wedge, Rotation
        from pytom_volume import rotate
        from pytom.basic.fourier import fft, ifft
        from pytom_volume import complexRealMult
        from pytom.basic.filter import rotateWeighting

        dim = 24
        cent = int(dim/2)
        v = vol(dim,dim,dim)
        v.setAll(0.0)
        v.setV(1.,cent,cent,cent)
        rot = Rotation(90.,0.,30.)
        wangle = 30.

        wedgeInfo = Wedge(wedgeAngles=wangle, cutoffRadius=0.0, tiltAxis='Y', smooth=3.0)
        tmp = wedgeInfo.apply(v)
        vfil1 = vol(v.sizeX(),v.sizeY(),v.sizeZ())
        rotate(tmp, vfil1, rot[0], rot[1], rot[2])

        wedgeInfoRot = Wedge(wedgeAngles=wangle, cutoffRadius=0.0,tiltAxis=rot,smooth = 3.0)
        vfil2 = wedgeInfoRot.apply(v)
        #vfil2.write('xxx2.em')

        wrot = wedgeInfo.returnWedgeVolume(v.sizeX(), v.sizeY(), v.sizeZ(), False, rot)
        fv = fft(v)
        fvol3 = complexRealMult(fv,wrot)
        vfil3 = ifft(fvol3)
        vfil3.shiftscale(0.0,1/float(v.sizeX()*v.sizeY()*v.sizeZ()))
        #vfil3.write('xxx3.em')

        w = wedgeInfo.returnWedgeVolume(v.sizeX(), v.sizeY(), v.sizeZ(), False)
        wrot = rotateWeighting( weighting=w, z1=rot[0], z2=rot[1], x=rot[2], mask=None,
                                isReducedComplex=None,returnReducedComplex=True)
        fvol4 = complexRealMult(fv,wrot)
        vfil4 = ifft(fvol4)
        vfil4.shiftscale(0.0,1/float(v.sizeX()*v.sizeY()*v.sizeZ()))
        #vfil4.write('xxx4.em')

        self.assertAlmostEqual(vfil1.getV(17,17,17),vfil2.getV(17,17,17),2,"Wedge Filter Rotation :(")
        self.assertAlmostEqual(vfil2.getV(17,17,17),vfil3.getV(17,17,17),2,"Wedge Filter Rotation :(")


    def test_correlation(self):
        """
        Why is this one here?
        """
        import pytom_volume
        from pytom.basic import correlation
        s = pytom_volume.vol(32,32,32)
        s2 = pytom_volume.vol(32,32,32)
        
        pytom_volume.initSphere(s,16,0,16,16,16,16)
        pytom_volume.initSphere(s2,16,0,16,16,16,16)
        
        xcf = correlation.xcf(s,s2)
        
        p = pytom_volume.peak(xcf)
        
        self.assertTrue(p[0] == 16 and p[1] == 16 and p[2] == 16, 
            msg='correlation test failed')


    def test_bandpassFilter(self):

        from pytom_volume import vol
        from pytom.basic.filter import bandpassFilter

        refNoise = vol(64,64,64)
        refNoise.setAll(0.)
        refNoise.setV(1., 32, 32, 32)
        # make sure that out of bounds frequency works
        fil1 = bandpassFilter(volume=refNoise, lowestFrequency=0, highestFrequency=30, bpf=None, smooth=1,
                             fourierOnly=False)
        fil2 = bandpassFilter(volume=refNoise, lowestFrequency=0, highestFrequency=16, bpf=None, smooth=1,
                             fourierOnly=False)
        #fil1[0].write('fil1.em')
        #fil2[0].write('fil2.em')
        filter = fil1[1]
        filter.getWeightVolume(True)

    
    def test_ones(self):
        """
	    output of filtering should be almost equal to input
        """
        import pytom_volume
        from pytom.basic.filter import filter
        from pytom_freqweight import weight
        
        vol1 = pytom_volume.vol(64,64,64)
        vol2 = pytom_volume.vol(64,64,64)
        vol1.setAll(1)
        vol1.setV(0,33,33,33)
			            
        vol2.setAll(1)
				             
        w = weight(vol2)
					              
        f = filter(vol1,w)
						       
        dif = vol1 - f[0]
        pytom_volume.abs(dif)
        p = pytom_volume.peak(dif)
        v = dif.getV(p[0],p[1],p[2])
        self.assertTrue( v < 0.0001, msg='ouput not equal to input')


    def test_ramp(self):
        """
        output of ramp filtering
        """
        from pytom_volume import vol, complexRealMult, peak
        from pytom.basic.filter import filter as fil
        from pytom.basic.filter import circleFilter,rampFilter,fourierFilterShift
        from pytom.basic.fourier import fft,ifft
        import os
        im = vol(64,64,1)
        
        im.setAll(0);
        im.setV(1,32,32,0);

        w = fourierFilterShift(rampFilter(64,64))
        
        weiproj = ifft( complexRealMult( fft( im), w) )/(im.sizeX()*
	          im.sizeY()*im.sizeZ())

        pos = peak(weiproj)
        val = weiproj.getV(pos[0],pos[1],pos[2])

        self.assertTrue( (pos[0] == 32) and (pos[1] == 32) and (pos[2] == 0), 
	    'rampFunction shifts pixels !')
        self.assertTrue( val==.5, 'original point too weak')

        # now test filter with bandpass
        w = fourierFilterShift(rampFilter(64,64))
        w.write('wei.em')

        weiproj = ifft( complexRealMult( fft( im), w) )/(im.sizeX()*
	          im.sizeY()*im.sizeZ())
        weiproj.write('weiproj.em')

        pos = peak(weiproj)
        val = weiproj.getV(pos[0],pos[1],pos[2])

        self.assertTrue( (pos[0] == 32) and (pos[1] == 32) and (pos[2] == 0),
	    'rampFunction shifts pixels !')
        self.assertTrue( val==.5, 'original point too weak')
        os.system('rm wei.em')
        os.system('rm weiproj.em')
        
    def test_complexDiv(self):
        from pytom_volume import vol, complexDiv,abs,peak
        from pytom.basic.fourier import fft,ifft
        
        v = vol(64,64,64)
        w = vol(64,64,33)
        sizeX = v.sizeX()
        sizeY = v.sizeY()
        sizeZ = v.sizeZ()
        
        v.setAll(0)
        v.setV(1,33,33,33)
        w.setAll(1)
        
        fv = fft(v);
        r = complexDiv(fv,w)
        result = ifft(r)       
        result.shiftscale(0.0,1/float(sizeX*sizeY*sizeZ))
        
        dif = v - result
        abs(dif)
        p = peak(dif)
        v = dif.getV(p[0],p[1],p[2])
        
        self.assertTrue( v < 0.00001, msg='')


    def test_profile(self):
        """test filter by 1-d profile"""
        from pytom.basic.filter import profile2FourierVol
        from pytom.basic.fourier import convolute, powerspectrum
        from pytom_volume import vol
        from math import sqrt
        
        #generate radial Fourier Filter
        profile = vol(16,1,1)
        for ii in range(0,16):
            profile.setV(ii, ii, 0, 0)
        kernel = profile2FourierVol( profile=profile, dim=None, reduced=False)
        
        # point 
        volume = vol(32,32,32)
        volume.setAll(0.)
        volume.setV(1,16,16,16)
        
        #make point spread function
        outvol = convolute(volume, kernel,True)
        ps = powerspectrum(outvol)
        
        sf = sqrt(32**3)
        for ii in range(0,16):
            d = abs(sqrt(ps.getV(ii+16,16,16))*sf - profile.getV(ii,0,0))
            self.assertTrue( d < 0.0001, 'ProfileTest broken :(')
        
             
    def runTest(self):
        self.test_wedgeFilter()
        self.test_wedgeRotation()
        self.test_bandpassFilter()
        self.test_ones()
        self.test_ramp()
        self.test_complexDiv()
        self.test_profile()
        
if __name__ == '__main__':
    unittest.main()


