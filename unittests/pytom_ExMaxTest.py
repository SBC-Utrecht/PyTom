import unittest
from pytom.basic.structures import ParticleList, Particle, Rotation, Shift, Wedge
from pytom_volume import vol, read, transformSpline
from random import seed, uniform
from pytom.alignment.alignmentFunctions import average

class pytom_ParticleListTest(unittest.TestCase):

    def setUp(self):
        """
        list with n ribosomes
        """
        from pytom.basic.structures import Reference, Mask, SampleInformation
        from pytom.basic.functions import initSphere
        from pytom_volume import gaussianNoise

        seed(0)
        self.npart=40
        self.angerr = 15 # initial angular assignment error
        self.noiselevel = 2
        #self.pl_name = 'testData/particleList.xml'
        self.templateName = 'testData/emd_1480.map.em_bin_4.em'
        self.average_name = 'testData/average.em'
        self.mask_name = 'testData/testMask.em'
        self.reference = Reference(referenceFile=self.templateName, generatedByParticleList=None)
        self.sample_info = SampleInformation(pixelSize=10.92, particleDiameter=250.)
        self.rot_ini = [0,0,0]
        self.trans_ini = [0,0,0]
        self.rot_width = [40.,40.,40.]
        self.trans_width = [2., 2., 2.]
        self.wedge = Wedge( wedgeAngles=[30.,30.], cutoffRadius=.0, tiltAxis='Y', smooth=2.0)
        self.template = read(self.templateName)
        maskvol = initSphere(sizeX=self.template.sizeX(), sizeY=self.template.sizeY(), sizeZ=self.template.sizeZ(),
                             radius=16, smooth=2, maxradius=20, cent=None, filename=self.mask_name)
        self.mask = Mask( filename=self.mask_name)
        self.voltrans = vol( self.template.sizeX(), self.template.sizeY(), self.template.sizeZ())
        noisevol = vol( self.template.sizeX(), self.template.sizeY(), self.template.sizeZ())
        self.pl_true = ParticleList()
        self.pl = ParticleList()
        for ipart in range(0, self.npart):
            rot = Rotation(uniform(self.rot_ini[0]-self.rot_width[0], self.rot_ini[0]+self.rot_width[0]),
                            uniform(self.rot_ini[1]-self.rot_width[1], self.rot_ini[1]+self.rot_width[1]),
                            uniform(self.rot_ini[2]-self.rot_width[2], self.rot_ini[2]+self.rot_width[2]) )
            trans = Shift(uniform(self.trans_ini[0]-self.trans_width[0], self.trans_ini[0]+self.trans_width[0]),
                            uniform(self.trans_ini[1]-self.trans_width[1], self.trans_ini[1]+self.trans_width[1]),
                            uniform(self.trans_ini[2]-self.trans_width[2], self.trans_ini[2]+self.trans_width[2]))
            transformSpline( self.template, self.voltrans, rot[0], rot[1], rot[2],
                    int(self.template.sizeX()/2),int(self.template.sizeY()/2),int(self.template.sizeY()/2),
                    0, 0, 0,
                    trans[0], trans[1], trans[2])
            fname = 'testData/particle_'+str(ipart)+'.em'
            cvol = self.wedge.apply( volume=self.voltrans)
            cvol.write(fname)
            p_true = Particle( filename=fname,rotation=rot, shift=trans, wedge=self.wedge)
            self.pl_true.append(particle=p_true)

            gaussianNoise(noisevol, 0., self.noiselevel)
            self.voltrans = self.voltrans + noisevol
            fname = 'testData/particleNoisy_'+str(ipart)+'.em'
            cvol = self.wedge.apply( volume=self.voltrans)
            cvol.write(fname)
            p = Particle( filename=fname,rotation=Rotation(uniform(rot[0]-self.angerr, rot[0]+self.angerr),
                                                           uniform(rot[1]-self.angerr, rot[1]+self.angerr),
                                                           uniform(rot[2]-self.angerr, rot[2]+self.angerr)),
                          shift=Shift(self.trans_ini[0],self.trans_ini[1],self.trans_ini[2]),
                          wedge=self.wedge)
            self.pl.append(particle=p)


    def average_Test(self):
        from pytom.basic.correlation import nxcc

        av_wei = average( particleList=self.pl_true, averageName=self.average_name, showProgressBar=False, verbose=False,
                createInfoVolumes=True, weighting=False, norm=False)
        av_wei = read(self.average_name)
        print("######### Averaging Test #########")
        av = read('testData/average-PreWedge.em')
        cc_wei = nxcc(volume=av_wei, template=self.template)
        print("CC(wedgeCorrected,orig)=", cc_wei)
        cc = nxcc(volume=av, template=self.template, mask=self.mask.getVolume())
        print("CC(unCorrected,orig)   =", cc)
        self.assertTrue(cc_wei > cc, "Fourier weighting results in lower CC than un-filtered FSC :(")
        p = self.pl_true[14]
        pvol = p.getTransformedVolume()
        cc = nxcc(volume=self.wedge.apply(volume=pvol, rotation=p.getRotation().invert()),
                  template=self.wedge.apply(volume=self.template, rotation=p.getRotation().invert()),
                  mask=self.mask.getVolume())
        print("CC(p[14],orig)   =", cc)


    def reference_Test(self):
        """
        test reference functions, in particular subtractAverage
        """
        from pytom.basic.structures import Reference
        from pytom.basic.correlation import nxcc

        ii = 19
        p = self.pl[ii]
        self.reference = Reference( referenceFile=self.average_name, generatedByParticleList=self.pl)
        [refVol,newWedgeSum] = self.reference.subtractParticle(particle=p, binning=1)

        print("######### Reference Subtraction Test #########")
        cc_wei = nxcc(refVol, self.template)
        print 'cc(correctedRef, template)=', cc_wei
        pvol = p.getTransformedVolume()
        cc1 = nxcc(self.wedge.apply(refVol), self.wedge.apply(pvol))
        print 'cc(correctedRef, particle)=', cc1
        cc2 = nxcc(self.wedge.apply(self.template), self.wedge.apply(pvol))
        print 'cc(template, particle)=', cc2
        self.assertTrue(cc2 > cc1, "Subtraction of particle from average results in mess :(")
        #refVol.write('testData/refVol.em')

    def exMax_Test(self):
        from pytom.alignment.ExMaxAlignment import ExMaxJob, sequentialStart
        from pytom.score.score import FLCFScore
        from pytom.angles.localSampling import LocalSampling
        from pytom.alignment.preprocessing import Preprocessing
        from pytom.basic.structures import Reference
        from pytom_volume import read

        self.reference = Reference(referenceFile=self.average_name, generatedByParticleList=self.pl)
        self.ejob = ExMaxJob(particleList=self.pl,destination='./testData/',reference=self.reference,
                 score=FLCFScore(), rotations=LocalSampling(shells=2,increment=3., z1Start=0., z2Start=0., xStart=0.),
                 mask=self.mask,symmetry=None,
                 numberRefinementRounds=3,numberIterations=3,preprocessing=Preprocessing(),
                 excludeThreshold=-1,binning=1,sampleInformation=self.sample_info, fscCriterion=0.5,
                 adaptiveResolution=True,adaptiveOffset=0.0,angleFactor=0.1)
        self.ejob.toXMLFile(filename='testData/exMaxJob.xml')

    def exMaxRun_Test(self):
        from pytom.alignment.ExMaxAlignment import sequentialStart
        sequentialStart(exMaxJob=self.ejob, verbose=False)


    def compareVolumes(self):
        #compare final average to template
        from pytom.basic.correlation import FSC
        av = read("testData/3.em")
        fsc = FSC(volume1=self.template, volume2=av, numberBands=av.sizeX()/2)
        print(fsc)



    def runTest(self):
        self.average_Test()
        self.reference_Test()
        self.exMax_Test()
        self.exMaxRun_Test()
        self.compareVolumes()



if __name__ == '__main__':
    unittest.main()
