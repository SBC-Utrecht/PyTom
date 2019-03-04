import unittest
epsilon = 0.001

class pytom_CorrelationTest(unittest.TestCase):
    
    def XCF_Test(self, a=0):
                
        import pytom_volume
        from pytom.basic import correlation
        from pytom.basic.functions import initSphere

        if a.__class__ == int: 
            a = pytom_volume.read('./testData/ribo.em')
            b = pytom_volume.read('./testData/ribo.em')

        mask = initSphere(sizeX=a.sizeX(), sizeY=a.sizeY(), sizeZ=a.sizeZ(),
	        radius=16,smooth=3,maxradius=a.sizeX()/2, 
		cent=[a.sizeX()/2, a.sizeY()/2, a.sizeZ()/2])
        
        xcfVol = correlation.xcf(volume=a, template=b, mask=mask)
        
        p = pytom_volume.peak(xcfVol)
        
        self.assertTrue( (p[0] == a.sizeX()/2) and (p[1] == a.sizeY()/2) 
	        and (p[2] == a.sizeZ()/2), 'XCF peak at wrong position')
        
        val = xcfVol.getV(p[0],p[1],p[2])
        self.assertTrue(152801232.0 + 152801232.0/5000000. >= val > 152801232.0 - 152801232.0/5000000.,
                        msg='Autocorrelation of XCF is wrong')
    
    def NXCC_Test(self,a=0):

        import pytom_volume
        from pytom.basic import correlation
        from pytom.basic.functions import initSphere
        
        if a.__class__ == int: 
            a = pytom_volume.read('./testData/ribo.em')
            b = pytom_volume.read('./testData/ribo.em')
        
        mask = initSphere(sizeX=a.sizeX(), sizeY=a.sizeY(), sizeZ=a.sizeZ(),
	        radius=16,smooth=3,maxradius=a.sizeX()/2, 
		cent=[a.sizeX()/2, a.sizeY()/2, a.sizeZ()/2])

        val = correlation.nxcc(volume=a, template=b, mask=mask)
        
        #print("NXCC="+str(val) )
        self.assertTrue(1.0 + epsilon >= val > 1.0 - epsilon, msg='Autocorrelation of NXCF is wrong')

    def NXCF_Test(self,a=0):

        import pytom_volume
        from pytom.basic import correlation
        from pytom.basic.functions import initSphere
        
        if a.__class__ == int: 
            a = pytom_volume.read('./testData/ribo.em')
            b = pytom_volume.read('./testData/ribo.em')
        
        mask = initSphere(sizeX=a.sizeX(), sizeY=a.sizeY(), sizeZ=a.sizeZ(),
	        radius=16,smooth=3,maxradius=a.sizeX()/2, 
		cent=[a.sizeX()/2, a.sizeY()/2, a.sizeZ()/2])
        xcfVol = correlation.nXcf(volume=a, template=b, mask=None)
        
        p = pytom_volume.peak(xcfVol)
        
        self.assertTrue( (p[0] == a.sizeX()/2) and (p[1] == a.sizeY()/2) 
	        and (p[2] == a.sizeZ()/2), 'NCF peak at wrong position')
	
        val = xcfVol.getV(p[0],p[1],p[2])
        #print("NXCF="+str(val) )
        self.assertTrue(1.0 + epsilon >= val > 1.0 - epsilon, msg='Autocorrelation of NXCF is wrong')
        
        from pytom.simulation.whiteNoise import add
        
        a = add(a)
        
        xcfVol = correlation.nXcf(a,b)
        
        p = pytom_volume.peak(xcfVol)
        
        assert p[0] == a.sizeX()/2 and p[1] == a.sizeY()/2 and p[2] == a.sizeZ()/2
        
        val = xcfVol.getV(p[0],p[1],p[2])
         
        assert 1 + epsilon >= val > -1 - epsilon

        
    def FLCF_Test(self,a=0):

        import pytom_volume

        # Test auto-correlation
        if a.__class__ == int: 
            a = pytom_volume.read('./testData/ribo.em')
            b = pytom_volume.read('./testData/ribo.em')
        
        m = pytom_volume.vol(a.sizeX(),a.sizeY(),a.sizeZ())
        pytom_volume.initSphere(m,a.sizeX()/2-1,0,0, 
	      a.sizeX()/2, a.sizeX()/2, a.sizeX()/2)
        
        b = b*m
        
        from pytom.basic.correlation import FLCF
        
        c = FLCF(a,b,m)
        
        p = pytom_volume.peak(c)
        v = c.getV(p[0],p[1],p[2])
        
        assert (1 + 0.01) >= v >= (1 - 0.01)
        
        # Rotate and test
        from pytom_volume import  rotateSpline
        rotateSpline(a, b, 30, 30, 30)
        c = FLCF(a,b,m)
        p = pytom_volume.peak(c)
        v = c.getV(p[0],p[1],p[2])
        
        assert 1  >= v >= -1 
        
        # Test fast local correlation function
        from pytom_volume import subvolume
        b = subvolume(a, 40,40,40, 32,32,32)
        c = FLCF(a,b)
        p = pytom_volume.peak(c)
        v = c.getV(p[0],p[1],p[2])
        
        assert (1 + 0.01) >= v >= (1 - 0.01)

    def SOC_Test(self,a=0):

        import pytom_volume;

        if a.__class__ == int: 
            a = pytom_volume.read('./testData/ribo.em')
            b = pytom_volume.read('./testData/ribo.em')
        
        from pytom.basic.correlation import soc
        
        c = soc(a,b)
        
        p = pytom_volume.peak(c)
        
        assert p[0] == a.sizeX()/2 and p[1] == a.sizeY()/2 and p[2] == a.sizeZ()/2
        
        v = c.getV(p[0],p[1],p[2])
        
        assert (1 + 0.001) >= v >= (1 - 0.001)
        
    
    def WXCC_Test(self,a=0):

        import pytom_volume

        if a.__class__ == int: 
            a = pytom_volume.read('./testData/ribo.em')
            b = pytom_volume.read('./testData/ribo.em')
        
        from pytom.basic.correlation import weightedXCC
        from pytom.score.score import RScore
        score = RScore()
        
        score.initAttributes()
         
        c = weightedXCC(a,b,10)
        
        assert (1 + epsilon) >= c >= 0.9
    
    def WXCF_Test(self,a=0):

        import pytom_volume;

        if a.__class__ == int: 
            a = pytom_volume.read('./testData/ribo.em')
            b = pytom_volume.read('./testData/ribo.em')
        
        from pytom.basic.correlation import weightedXCF
        from pytom.score.score import RScore
        score = RScore()
        
        score.initAttributes(10)
         
        cc = weightedXCF(a,b,10)
        
        p = pytom_volume.peak(cc)
        c = cc.getV(p[0],p[1],p[2])
        assert (1 + epsilon) >= c >= (1 - epsilon)

    
    def FSC_Test(self):
 
        import pytom_volume
        from pytom.basic.correlation import FSC
        
        a = pytom_volume.vol(20,20,20)
        pytom_volume.gaussianNoise(a,0,1)
        
        b = pytom_volume.vol(20,20,20)
        b.copyVolume(a)
        
        c = pytom_volume.vol(20,20,20)
        pytom_volume.gaussianNoise(c,0,1)
        
        fsc = FSC(a,b,10)
        
        for i in range(len(fsc)-1):
            assert abs(1 - fsc[i]) <= epsilon
        
        fsc = FSC(a,c,10);
        
        for i in range(len(fsc)-1):
            assert  abs(abs(fsc[i]) - 1) <= 1
        
    def SubPixelPeak_T(self,volume,coordinates,resultCoordinates,resultValue,distance):    
        from pytom_volume import vol
        from pytom.basic.correlation import subPixelPeak
        from pytom.tools.maths import euclidianDistance
         
        r = subPixelPeak(volume, coordinates)

        assert r[0] >= resultValue
        
        assert euclidianDistance(r[1],resultCoordinates) <= distance
        
    def SubPixelPeak_Test(self):
        """
        
        """
        from pytom_volume import vol,transformSpline
        import pytom_volume
        
        v = vol(32,32,32)
        v.setAll(0)
        v.setV(1,12,12,12)
        
        self.SubPixelPeak_T(v,[12,12,12],[12,12,12],1,0)
        
        v.setAll(0)
        v.setV(1,12,13,14)
        
        self.SubPixelPeak_T(v,[12,13,14],[12,13,14],1,0)
        
        v.setAll(0)
        v.setV(1,8,8,8)
        v2 = vol(32,32,32)
        transformSpline(v,v2,0,0,0,0,0,0,0.5,0.5,0.5,0,0,0)

        self.SubPixelPeak_T(v2,[8,8,8],[8.4,8.4,8.4],0.125,0.2)
        
        a = pytom_volume.read('./testData/ribo.em')
        b = pytom_volume.read('./testData/ribo.em')
        
        from pytom.basic import correlation
        
        sizeX = a.sizeX()
        sizeY = a.sizeY()
        sizeZ = a.sizeZ()
        
        xcfVol = correlation.nXcf(a,b)
        
        self.SubPixelPeak_T(xcfVol,[sizeX//2,sizeY//2,sizeZ//2],[sizeX//2,sizeY//2,sizeZ//2],1,0.0)
        
        c = vol(sizeX,sizeY,sizeZ)
        transformSpline(b,c,0,0,0,0,0,0,2.2,5.6,4.8,0,0,0)
        
        xcfVol = correlation.nXcf(c,a)
        
        self.SubPixelPeak_T(xcfVol,[sizeX//2+2,sizeY//2+5,sizeZ//2+4],[sizeX//2+2.2,sizeY//2+5.6,sizeZ//2+4.8],0.99,epsilon)
        
    def runTest(self):
        self.XCF_Test()
        self.NXCC_Test()
        
if __name__ == '__main__':
    unittest.main()

