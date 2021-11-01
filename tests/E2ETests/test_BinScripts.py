import os
import sys
from pytom.agnostic.io import read, write
import numpy as np
import unittest
from numpy import sqrt
from os import system as exe
from os.path import join as merge
from pytom.basic.structures import Particle, ParticleList
from os.path import exists

class pytom_MyFunctionTest(unittest.TestCase):

    def setUp(self):
        from helper_functions import create_RandomParticleList, installdir

        self.installdir = installdir
        self.reffile = f'../testData/ribo.em'
        self.pdir = f'./testparticles'


        # set parameters for ACL
        self.settings = {}
        self.settings["frequency"] = 40
        self.settings["outputDirectory"] = 'outputCCC/'
        if not os.path.exists(self.settings["outputDirectory"]):
            os.mkdir(self.settings["outputDirectory"])

        self.pl_filename = os.path.join(self.settings['outputDirectory'], 'pl.xml')

        self.pl = create_RandomParticleList(reffile=self.reffile, pl_filename=self.pl_filename,
                                            pdir=self.pdir, nparticles=10)
        #self.settings["fixed_frequency"] = True
        #self.settings["offset"] = None
        #self.settings["mask"] = options.mask

        self.settings["mask"] = f'../testData/mask_45.em'
        self.folder = 'UnitTestFolder'
        self.create_RandomParticleList = create_RandomParticleList
        self.outdir = self.settings["outputDirectory"]

    def addProjectDirToParticleList(self):
        from pytom.basic.structures import ParticleList
        fname = 'dummy.xml'
        part = 'Structures/particle_0.mrc'
        folder = self.outdir
        pl = ParticleList()
        pl.append(Particle(part))
        pl.toXMLFile(fname)

        os.system(f'addProjectDirToParticleList.py -p {fname} -d {folder} -o {fname}')

        pl2 = ParticleList()
        pl2.fromXMLFile(fname)

        assert pl2[0].getFilename() == os.path.join(folder, part)

    def average(self):
        func = self.generate_cmd(sys._getframe().f_code.co_name)

        outname = 'average.mrc'
        self.pl.toXMLFile('average.xml')
        cmd  = f'average.py -p average.xml -a {outname} -c 1'

        self.check_cmd(cmd, func, outname)

    def bandpassFilterVolume(self):
        outname = 'dummy.mrc'
        os.system(f'bandpassFilterVolume.py -v {self.reffile} -o {outname} -l 10 -h 25 -s 3')
        self.assertTrue(os.path.exists(outname))
        os.remove(outname)

    def cancel_batch(self):
        '''Cancels all slurm jobs between first number and last number (including last)'''
        os.system('cancel_batch 10000 100001')

    def combineParticleLists(self):

        outdir = 'rest'
        self.create_folder(outdir)
        outname0 = os.path.join(outdir, 'dummy0.xml')
        outname1 = os.path.join(outdir, 'dummy1.xml')
        outfile = merge(outdir, 'combined.xml')
        pl = self.pl

        pl.toXMLFile(outname0)
        pl.toXMLFile(outname1)

        func = self.generate_cmd(sys._getframe().f_code.co_name)

        cmd = f'{func} -d {outdir} -o {outfile}'

        self.check_cmd(cmd, func, outfile)

    def check_cmd(self, cmd, func, outfile):
        print(cmd)
        exe(cmd)
        self.assertTrue(exists(outfile), f'{func} failed')
        try:
            os.remove(outfile)
        except:
            os.system(f'rm -rf {outfile}')

    def convert(self):
        folder = 'testparticles'
        func = self.generate_cmd(sys._getframe().f_code.co_name)
        cmd = f"{func} -d {folder} -t {self.outdir} -o mrc "
        self.check_cmd(cmd, func, merge(self.outdir, 'particle_0.mrc'))

    def create3DEllipse(self):
        outfile = merge(self.outdir, 'ellipse.mrc')
        func = self.generate_cmd(sys._getframe().f_code.co_name)
        cmd = f"{func} -e {128} -m {40} -n {30} -l {20} -o {outfile} -s {5} -c {3}"
        self.check_cmd(cmd, func, outfile)

    def createFolderStructure(self):
        outfolder = merge(self.outdir, 'guiFolderStructure')
        func = self.generate_cmd(sys._getframe().f_code.co_name)
        cmd = f"{func} -d {outfolder}"
        self.check_cmd(cmd, func, outfolder)

    def createParticleListFromDir(self):
        outfile = merge(self.outdir, 'createPL.xml')
        func = self.generate_cmd(sys._getframe().f_code.co_name)
        prefix = 'particle_'
        cmd = f"{func} -d {self.outdir} -p {prefix} -o {outfile} -w 10"

        self.check_cmd(cmd, func, outfile)

    def cropSubvolumes(self):
        outfile = merge(self.outdir, 'cropSubPL.xml')
        func = self.generate_cmd(sys._getframe().f_code.co_name)
        inname = merge(self.outdir, 'start.xml')
        self.pl.toXMLFile(inname)
        cmd = f"{func} --particleList {inname} --output particle_ --center 50,50,50 --cubesize 20 --outParticleList {outfile}"

        self.check_cmd(cmd, func, outfile)

    def diff(self):
        os.system(f'diff.py {self.pl[0].getFilename()} {self.pl[1].getFilename()}')

    def extractClassesFromParticleList(self):
        plname = merge(self.outdir, 'class.xml')

        pl = self.pl.copy()
        for n, particle in enumerate(pl):
            particle.setClass(f'{n}')

        pl.toXMLFile(plname)
        outfile = merge(self.outdir, 'sel.xml')
        func = self.generate_cmd(sys._getframe().f_code.co_name)
        cmd = f"{func} -p {plname}  -c 1,2 -o {outfile}"

        self.check_cmd(cmd, func, plname)

    def extractClassXML(self):
        plname = merge(self.outdir, 'class2.xml')

        pl = self.pl.copy()
        for n, particle in enumerate(pl):
            particle.setClass(f'{n}')

        pl.toXMLFile(plname)
        outfile = merge(self.outdir, 'sel.xml')
        func = self.generate_cmd(sys._getframe().f_code.co_name)
        cmd = f"{func} -p {plname}  -c 1,2 -o {outfile}"

        self.check_cmd(cmd, func, plname)

        pl = ParticleList()
        pl.fromXMLFile(outfile)
        for p in pl:
            self.assertTrue( p.getClass() in ("1", '2'), 'wrong class')

    def extractProjectDirFromParticleList(self):
        plname = merge(self.outdir, 'dd.xml')

        pl = self.pl.copy()
        for n, particle in enumerate(pl):
            particle.getInfoGUI().setProjectDir(f'Project_{0}')

        pl.toXMLFile(plname)
        outfolder = merge(self.outdir, 'particleList')
        self.create_folder(outfolder)
        func = self.generate_cmd(sys._getframe().f_code.co_name)
        cmd = f"{func} -p {plname}  -d {outfolder}"

        self.check_cmd(cmd, func, plname)

    def extractTomoNameFromXML(self):
        plname = merge(self.outdir, 'dd.xml')

        pl = self.pl.copy()
        for n, particle in enumerate(pl):
            particle.getPickPosition().setOriginFilename(f'tomogram_{n:03d}_WBP.mrc')

        pl.toXMLFile(plname)

        outfolder = merge(self.outdir, 'particleList')
        self.create_folder(outfolder)
        func = self.generate_cmd(sys._getframe().f_code.co_name)
        cmd = f"{func} -p {plname}  -t {outfolder}"

        self.check_cmd(cmd, func, plname)

    def filter(self):
        outname = merge(self.settings['outputDirectory'], 'dance.mrc')
        func = self.generate_cmd(sys._getframe().f_code.co_name)
        cmd = f'{func} -f {self.reffile} -t {outname} -l 10 --highestFrequency 25 -s 3'
        os.system(cmd)
        self.assertTrue(os.path.exists(outname), f'{func} failed')

    def flip_coordinate_file(self):
        coordinateFile = '../testData/coords_tomogram_000_WBP.txt'
        outCoordinateFile = 'coords_tomogram_000_WBP_flipped.txt'
        func = self.generate_cmd(sys._getframe().f_code.co_name)
        cmd = f"{func} -f {coordinateFile} -d ../testData -t {self.settings['outputDirectory']} -s 464 -p flipped"
        os.system(cmd)
        self.assertTrue(os.path.exists(f"{self.settings['outputDirectory']}/{outCoordinateFile}"), f'{func} failed')

    def fsc(self):
        func = self.generate_cmd(sys._getframe().f_code.co_name)
        fsc = 0.143
        outfolder = self.outdir
        pixelsize=10
        if not os.path.exists(outfolder): os.mkdir(outfolder)
        cmd = f'{func} --v1 {self.reffile} --v2 {self.reffile} --fsc {fsc} --pixelsize {pixelsize} --combinedResolution --outputFolder {outfolder} --xml'
        self.check_cmd(cmd, func, os.path.join(outfolder, 'FSCOrig.dat'))

    def gen_mask(self):
        outname = 'dummy.mrc'
        os.system(f'gen_mask.py -f {self.reffile} -o {outname} -n {1} -s {3} -c {4}')
        self.assertTrue(os.path.exists(outname))
        os.remove(outname)

    def lenPL(self):
        fname = merge(self.outdir, 'len.xml')
        self.pl.toXMLFile(fname)
        self.assertTrue('10' == os.popen(f'lenPL.py -f {fname}').read()[:-1].split()[-1], 'lenPL.py is not functional')

    def create_folder(self, fname):
        if not os.path.exists(fname): os.mkdir(fname)

    def mirrorVolume(self):
        outname = merge(self.outdir, 'dummy.em')
        os.system(f'mirrorVolume.py -v {self.reffile} -o {outname}')
        data = read(outname)
        ref = read(self.reffile)

        diff = np.abs(data[::-1,:,:] - ref)
        self.assertTrue(diff.sum() < 1E-5)

    def mrcs2mrc(self):
        from pytom.agnostic.io import read, write
        import numpy as np

        fname = 'dummy.mrcs'
        size =11
        vol = np.zeros((size,size,size),dtype=np.float32)
        for i in range(size):
            vol[i,i,i] = i+1

        write(fname, vol)

        os.system(f'mrcs2mrc.py -f {fname} -t {self.outdir} -p sorted_ -i 12')

        for i in range(size):
            self.assertTrue(os.path.exists(f'{self.outdir}/sorted_{i:02d}.mrc'), 'file not generated')
            data = read(f'{self.outdir}/sorted_{i:02d}.mrc').squeeze()
            
            self.assertTrue(np.abs(data[i,i] - (i+1))  < 1E-5, 'Max value is off')
            self.assertTrue(np.abs(data).sum() > i+1-1E-5, f'{np.abs(data).sum()} ')

    def setWedgeToParticleList(self):
        angle = 45
        outfile = merge(self.outdir, 'addedWedge.xml')
        plname = merge(self.outdir, 'temp.xml')
        self.pl.toXMLFile(plname)
        func = self.generate_cmd(sys._getframe().f_code.co_name)
        cmd = f"{func} -p {plname} -o {outfile} -w {angle}"
        os.system(cmd)

        pl = ParticleList()
        pl.fromXMLFile(outfile)
        for particle in pl:
            self.assertTrue(np.abs(particle.getWedge().getWedgeAngle() - angle) < 1E-5,
                            f'wedge angle is not set to {angle}')

    def symmetrize(self):
        plname = merge(self.outdir, 'temp.xml')
        self.pl.toXMLFile(plname)
        volume = 'testparticles/particle_0.em'
        outfile = merge(self.outdir, 'cropSubPL.xml')
        func = self.generate_cmd(sys._getframe().f_code.co_name)
        cmd = f"{func} -p {plname} -v {volume} -r {outfile} -s {3}"

        self.check_cmd(cmd, func, outfile)


    # Already in Micrograph modeller unittest
    def template_generation(self):
        pass

    # Need to find solution to test plotting, maybe flag that not show but saves fig

    def mrcshow(self):
        pass

    def plotAngularDistribution(self):
        pass

    def plotGaussianFit(self):
        pass



    # here we need a good working example. It is duplicate with FullPipeline, so kept it to pass.

    def ctfCorrImod(self):
        pass

    def reconstructTomogramGPU(self):
        pass

    def reconstructTomogram(self):
        pass

    def calcAngThickSpecimen(self):
        pass

    def reconstructTomogramWithoutAlignment(self):
        pass

    def localizationJob(self):
        folder = ''
        dest = self.settings['outputDirectory']
        self.create_folder(dest)

        cmd = self.generate_cmd(sys._getframe().f_code.co_name)

        cmd = f'cd {self.folder}; mpiexec -np 16 {cmd} -j {jobname} '
        cmd += f'--volume {self.reffile} '
        cmd += f'--reference {self.reffile} '
        cmd += f'--mask {self.mask} '
        cmd += f'--wedge1 {30} '
        cmd += f'--wedge2 {30} '
        cmd += f'--angles angles_50_100.em '
        cmd += f'--destination {dest} '
        cmd += f'-b 40 '
        cmd += f'--splitX 1 --splitY 1 --splitZ 1'

        os.system(cmd)

    def localization(self):

        jobname = 'testjob.xml'

        cmd = self.generate_cmd(sys._getframe().f_code.co_name)

        cmd = f'cd {self.folder}; mpiexec -np 16 {cmd} -x 4 -y 4 -z 1 -j {jobname}'

        os.system(cmd)

    def extractCandidates(self):
        pass

    def templateMatchingCandidateExtractionSingleGPU(self):
        pass

    def cutParticlesFromVol(self):
        pass

    def deconv(self):
        pass

    def reconstructWB(self):
        pass

    def updateParticleList(self):
        pass

    def reconstructWB2(self):
        pass

    def alignJob(self):
        pass

    def align(self):
        pass

    def GLocalJob(self):
        pass

    def CPCAJob(self):
        pass

    def mcoACJob(self):
        pass

    def mcoAC(self):
        pass

    def mcoEXMXJob(self):
        pass

    def mcoEXMX(self):
        pass



    def generate_cmd(self, name):

        return(name + '.py ')

    def runTest(self):
        self.addProjectDirToParticleList()
        self.average()
        self.bandpassFilterVolume()
        self.cancel_batch()
        self.combineParticleLists()
        self.convert()
        self.create3DEllipse()
        self.createFolderStructure()
        self.createParticleListFromDir()
        self.cropSubvolumes()
        self.diff()
        self.extractClassesFromParticleList()
        self.extractClassXML()
        self.extractProjectDirFromParticleList()
        self.extractTomoNameFromXML()
        self.filter()
        self.flip_coordinate_file()
        self.fsc()
        self.gen_mask()
        self.lenPL()
        self.mirrorVolume()
        self.mrcs2mrc()
        self.setWedgeToParticleList()
        # self.symmetrize()

if __name__ == '__main__':
    unittest.main()
