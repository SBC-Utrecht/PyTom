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

    def test_addProjectDirToParticleList(self):
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
                    
        assert pl2[0].getInfoGUI().getProjectDir() == folder
        assert pl2[0].getFilename() == part

    def test_average(self):
        func = 'average.py'

        outname = merge(self.outdir, 'average.mrc')
        self.pl.toXMLFile('average.xml')
        cmd  = f'average.py -p average.xml -a {outname} -c 1'

        self.check_cmd(cmd, func, outname)

    def test_bandpassFilterVolume(self):
        outname = merge(self.outdir, 'dummy.mrc')
        os.system(f'bandpassFilterVolume.py -v {self.reffile} -o {outname} -l 10 -h 25 -s 3')
        self.assertTrue(os.path.exists(outname))
        os.remove(outname)

    def test_cancel_batch(self):
        '''Cancels all slurm jobs between first number and last number (including last)'''
        os.system('cancel_batch.py 10000 100001')

    def test_combineParticleLists(self):

        outdir = 'rest'
        self.create_folder(outdir)
        outname0 = os.path.join(outdir, 'dummy0.xml')
        outname1 = os.path.join(outdir, 'dummy1.xml')
        outfile = merge(outdir, 'combined.xml')
        pl = self.pl

        pl.toXMLFile(outname0)
        pl.toXMLFile(outname1)

        func = 'combineParticleLists.py'

        cmd = f'{func} -d {outdir} -o {outfile}'

        self.check_cmd(cmd, func, outfile)

    def check_cmd(self, cmd, func, outfile):
        exe(cmd)
        self.assertTrue(exists(outfile), f'{func} failed')
        try:
            os.remove(outfile)
        except:
            os.system(f'rm -rf {outfile}')

    def test_convert(self):
        folder = 'testparticles'
        func = 'convert.py'
        cmd = f"{func} -d {folder} -t {self.outdir} -o mrc "
        self.check_cmd(cmd, func, merge(self.outdir, 'particle_0.mrc'))

    def test_create_mask(self):
        outfile = merge(self.outdir, 'ellipse.mrc')
        func = 'create_mask.py'
        cmd = f"{func} --boxSize {128} --radius {40} --minor1 {30} --minor2 {20} -o {outfile} --sigma {5} --cutoff {3}"
        self.check_cmd(cmd, func, outfile)

    def test_createFolderStructure(self):
        outfolder = merge(self.outdir, 'guiFolderStructure')
        func = 'createFolderStructure.py'
        cmd = f"{func} -d {outfolder}"
        self.check_cmd(cmd, func, outfolder)

    def test_createParticleListFromDir(self):
        outfile = merge(self.outdir, 'createPL.xml')
        func = 'createParticleListFromDir.py'
        prefix = 'particle_'
        cmd = f"{func} -d {self.outdir} -p {prefix} -o {outfile} -w 10"

        self.check_cmd(cmd, func, outfile)

    def test_cropSubvolumes(self):
        outfile = merge(self.outdir, 'cropSubPL.xml')
        func = 'cropSubvolumes.py'
        inname = merge(self.outdir, 'start.xml')
        self.pl.toXMLFile(inname)
        cmd = f"{func} --particleList {inname} --output particle_ --center 50,50,50 --cubesize 20 --outParticleList {outfile}"

        self.check_cmd(cmd, func, outfile)

    def test_diff(self):
        os.system(f'diff.py {self.pl[0].getFilename()} {self.pl[1].getFilename()}')

    def test_extractClassesFromParticleList(self):
        plname = merge(self.outdir, 'class.xml')

        pl = self.pl.copy()
        for n, particle in enumerate(pl):
            particle.setClass(f'{n}')

        pl.toXMLFile(plname)
        outfile = merge(self.outdir, 'sel.xml')
        func = 'extractClassesFromParticleList.py'
        cmd = f"{func} -p {plname}  -c 1,2 -o {outfile}"

        self.check_cmd(cmd, func, plname)

    def test_extractClassXML(self):
        plname = merge(self.outdir, 'class2.xml')

        pl = self.pl.copy()
        for n, particle in enumerate(pl):
            particle.setClass(f'{n}')

        pl.toXMLFile(plname)
        outfile = merge(self.outdir, 'sel.xml')
        func = 'extractClassXML.py'
        cmd = f"{func} -p {plname}  -c 1,2 -o {outfile}"
        self.check_cmd(cmd, func, plname)

        pl = ParticleList()
        #Mirror what the bin script does
        real_out = outfile.replace('.xml','_deselected.xml')
        pl.fromXMLFile(real_out)
        for p in pl:
            self.assertTrue( p.getClass() in ("1", '2'), 'wrong class')

    def test_extractProjectDirFromParticleList(self):
        plname = merge(self.outdir, 'dd.xml')

        pl = self.pl.copy()
        for n, particle in enumerate(pl):
            particle.getInfoGUI().setProjectDir(f'Project_{0}')

        pl.toXMLFile(plname)
        outfolder = merge(self.outdir, 'particleList')
        self.create_folder(outfolder)
        func = 'extractProjectDirFromParticleList.py'
        cmd = f"{func} -p {plname}  -d {outfolder}"

        self.check_cmd(cmd, func, plname)

    def test_extractTomoNameFromXML(self):
        plname = merge(self.outdir, 'dd.xml')

        pl = self.pl.copy()
        for n, particle in enumerate(pl):
            particle.getPickPosition().setOriginFilename(f'tomogram_{n:03d}_WBP.mrc')

        pl.toXMLFile(plname)

        outfolder = merge(self.outdir, 'particleList')
        self.create_folder(outfolder)
        func = 'extractTomoNameFromXML.py'
        cmd = f"{func} -p {plname}  -t {outfolder}"

        self.check_cmd(cmd, func, plname)

    def test_filter(self):
        outname = merge(self.settings['outputDirectory'], 'dance.mrc')
        func = 'filter.py'
        cmd = f'{func} -f {self.reffile} -t {outname} -l 10 --highestFrequency 25 -s 3'
        os.system(cmd)
        self.assertTrue(os.path.exists(outname), f'{func} failed')

    def test_flip_coordinate_file(self):
        coordinateFile = '../testData/coords_tomogram_000_WBP.txt'
        outCoordinateFile = 'coords_tomogram_000_WBP_flipped.txt'
        func = 'flip_coordinate_file.py'
        cmd = f"{func} -f {coordinateFile} -d ../testData -t {self.settings['outputDirectory']} -s 464 -p flipped"
        os.system(cmd)
        self.assertTrue(os.path.exists(f"{self.settings['outputDirectory']}/{outCoordinateFile}"), f'{func} failed')

    def test_fsc(self):
        func = 'fsc.py'
        fsc = 0.143
        outfolder = self.outdir
        pixelsize=10
        if not os.path.exists(outfolder): os.mkdir(outfolder)
        cmd = f'{func} --v1 {self.reffile} --v2 {self.reffile} --fsc {fsc} --pixelsize {pixelsize} --combinedResolution --outputFolder {outfolder} --xml'
        self.check_cmd(cmd, func, os.path.join(outfolder, 'FSCOrig.dat'))

    def test_gen_mask(self):
        outname = merge(self.outdir, 'dummy.mrc')
        os.system(f'gen_mask.py -f {self.reffile} -o {outname} -n {1} -s {3} -c {4}')
        self.assertTrue(os.path.exists(outname))
        os.remove(outname)

    def test_lenPL(self):
        fname = merge(self.outdir, 'len.xml')
        self.pl.toXMLFile(fname)
        self.assertTrue('10' == os.popen(f'lenPL.py -f {fname}').read()[:-1].split()[-1], 'lenPL.py is not functional')

    def create_folder(self, fname):
        if not os.path.exists(fname): os.mkdir(fname)

    def test_mirrorVolume(self):
        outname = merge(self.outdir, 'dummy.em')
        os.system(f'mirrorVolume.py -v {self.reffile} -o {outname}')
        data = read(outname)
        ref = read(self.reffile)

        diff = np.abs(data[::-1,:,:] - ref)
        self.assertTrue(diff.sum() < 1E-5)

    def test_mrcs2mrc(self):
        from pytom.agnostic.io import read, write
        import numpy as np

        fname = merge(self.outdir, 'dummy.mrcs')
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

    def test_setWedgeToParticleList(self):
        angle = 45
        outfile = merge(self.outdir, 'addedWedge.xml')
        plname = merge(self.outdir, 'temp.xml')
        self.pl.toXMLFile(plname)
        func = 'setWedgeToParticleList.py'
        cmd = f"{func} -p {plname} -o {outfile} -w {angle}"
        os.system(cmd)

        pl = ParticleList()
        pl.fromXMLFile(outfile)
        for particle in pl:
            self.assertTrue(np.abs(particle.getWedge().getWedgeAngle() - angle) < 1E-5,
                            f'wedge angle is not set to {angle}')

    def test_symmetrize(self):
        plname = merge(self.outdir, 'temp.xml')
        self.pl.toXMLFile(plname)
        volume = 'testparticles/particle_0.em'
        outfile = merge(self.outdir, 'cropSubPL.xml')
        func = "symmetrize.py"
        cmd = f"{func} -p {plname} -v {volume} -r {outfile} -s {3}"
        self.check_cmd(cmd, func, outfile)


    # Already in Micrograph modeller unittest
    def test_create_template(self):
        pass

    # Need to find solution to test plotting, maybe flag that not show but saves fig

    def test_mrcshow(self):
        pass

    def test_plotAngularDistribution(self):
        pass

    def test_plotGaussianFit(self):
        pass



    # here we need a good working example. It is duplicate with FullPipeline, so kept it to pass.

    def test_ctfCorrImod(self):
        pass

    def test_reconstructTomogramGPU(self):
        pass

    def test_reconstructTomogram(self):
        pass

    def test_calcAngThickSpecimen(self):
        pass

    def test_reconstructTomogramWithoutAlignment(self):
        pass

    def test_localizationJob(self):
        folder = ''
        jobname = 'testlocaljob.xml'
        dest = self.settings['outputDirectory']
        self.create_folder(dest)
        mpi_cores = 4
        cmd = 'localizatonJob.py'
        mask = self.settings["mask"]

        cmd = f'cd {self.folder}; mpiexec -np {mpi_cores} {cmd} -j {jobname} '
        cmd += f'--volume {self.reffile} '
        cmd += f'--reference {self.reffile} '
        cmd += f'--mask {mask} '
        cmd += f'--wedge1 {30} '
        cmd += f'--wedge2 {30} '
        cmd += f'--angles angles_50_100.em '
        cmd += f'--destination {dest} '
        cmd += f'-b 40 '
        cmd += f'--splitX 1 --splitY 1 --splitZ 1'

        os.system(cmd)

    def test_localization(self):
        # TODO: deal with this test
        raise unittest.SkipTest('testjob.xml has disapeared, should be recovered')
        jobname = 'testjob.xml'
        mpi_procs = 4
        cmd = 'localization.py'

        cmd = f'cd {self.folder}; mpiexec -np {mpi_procs} {cmd} -x 4 -y 4 -z 1 -j {jobname}'

        os.system(cmd)

    def test_extractCandidates(self):
        pass

    def test_templateMatchingCandidateExtractionSingleGPU(self):
        pass

    def test_cutParticlesFromVol(self):
        pass

    def test_deconv(self):
        pass

    def test_reconstructWB(self):
        pass

    def test_updateParticleList(self):
        pass

    def test_reconstructWB2(self):
        pass

    def test_alignJob(self):
        pass

    def test_align(self):
        pass

    def test_GLocalJob(self):
        pass

    def test_CPCAJob(self):
        pass

    def test_mcoACJob(self):
        pass

    def test_mcoAC(self):
        pass

    def test_mcoEXMXJob(self):
        pass

    def test_mcoEXMX(self):
        pass


if __name__ == '__main__':
    unittest.main()
