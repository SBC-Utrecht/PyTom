import unittest
import os
from pytom.gpu.initialize import device
import sys

print(device)

class pytom_MyFunctionTest(unittest.TestCase):

    def setUp(self):
        from helper_functions import create_RandomParticleList
        import os

        pythonversion = f'python{sys.version_info.major}.{sys.version_info.minor}'
        if pythonversion == 'python3.7': pythonversion += 'm'


        self.projectname = f'{os.getcwd()}/FullPipeline'
        self.refDataDir = f'{os.getcwd()}/../testData'
        self.tomoname = f'{self.projectname}/03_Tomographic_Reconstruction/tomogram_000'
        self.metafile = f'{self.tomoname}/sorted/mixedCTEM_tomo3.meta'
        self.markerFile = '../testData/markerfile.txt'

        self.orignal_data = 'ftp://ftp.ebi.ac.uk/empiar/world_availability/10064/data/mixedCTEM_tomo3.mrc'
        self.mrcs = 'mixedCTEM_tomo3.mrc'
        self.pytomDir = os.path.dirname(os.path.dirname(os.popen('which pytom').readlines()[0]))
        if 'miniconda' in self.pytomDir:
            self.pytomDir2 = os.path.join(self.pytomDir, f'lib/{pythonversion}/site-packages/pytom')
        else:
            self.pytomDir2 = self.pytomDir

        self.alignment_sorted_results = f'../testData/alignmentResultsFullRange.txt'
        self.tmfolder =f'{self.projectname}/04_Particle_Picking/Template_Matching/cross_correlation/tomogram_000_WBP'
        self.plFilename = f'{self.projectname}/04_Particle_Picking/Picked_Particles/particleList_TM_tomogram_000_WBP_deselected.xml'
        self.plFilenameReduced = f'{self.projectname}/04_Particle_Picking/Picked_Particles/particleList_TM_tomogram_000_WBP_GLocal_reduced_best.xml'

        self.glocaldir = f'{self.projectname}/05_Subtomogram_Analysis/Alignment/GLocal/'

        self.singleGpuID = 1
        self.ctfFileNameAngle = f'{self.refDataDir}/angles.tlt'
        self.ctfFileNameDefocus = f'{self.refDataDir}/defocusValuesGijs.defocus'

        self.weightingTypeRecon = -1
        self.weightingTypeSubtomoRecon = -1
        self.weightingTypeAlignment = 0

        self.referenceMarker = 4
        self.referenceTiltIndex = 30
        self.expectedRotationAngle = 0
        self.firstIndexFull = 0
        self.lastIndexFull = 58
        self.firstIndexReduced = 20
        self.lastIndexReduced = 40
        self.startAngleFull = -60
        self.endAngleFull = 56
        self.startAngleReduced = -20
        self.endAngleReduced = 20

        self.numcores = 15
        self.IMODTiltAxis = 180

        self.pixelSize = 2.62
        self.binningGLocalFull = 2
        self.binningGLocalReduced = 1
        self.particleDiameter = 300
        self.dont = True
        global device

        device= 'gpu:0'

    def cleanUp(self):
        """
        check that files are written and remove them
        """
        from helper_functions import cleanUp_RandomParticleList
        import os
        #os.system(f'rm -rf {self.projectname}')

        # cleanUp_RandomParticleList(pl_filename=self.pl_filename, pdir=self.pdir)

    def remove_file(self, filename):
        """
        assert that file exists end remove it
        """
        from os import remove
        from os import path

        filecheck = path.exists(filename)
        self.assertTrue(filecheck, msg="file " + filename + " does not exist")
        if filecheck:
            remove(filename)

    def test_00_ProjectFolderGeneration(self):
        import os
        from pytom.gui.guiFunctions import create_project_filestructure
        create_project_filestructure(self.projectname)
        self.assertTrue( os.path.exists(self.projectname), msg="folder "+self.projectname+" does not exist")

        os.system(f'cp -rf {self.projectname}/03_Tomographic_Reconstruction/.tomoname {self.projectname}/03_Tomographic_Reconstruction/tomogram_000')
        pass

    def test_01_dataDownload(self):
        import os

        if not os.path.exists(self.mrcs):
            os.system(f'wget {self.orignal_data}')

        self.assertTrue(os.path.exists(self.mrcs), msg=f'Download of {self.mrcs} failed')
        pass

    def test_02_dataExtraction(self):
        import os
        if self.dont: return
        os.system(f'mrcs2mrc.py -f {self.mrcs} -t {self.tomoname}/sorted -p sorted -i 2 -s {self.startAngleFull} -e {self.endAngleFull} -m ')
        os.system(f'cp {self.markerFile} {self.tomoname}/sorted')

        self.assertTrue(os.path.exists(f'{self.tomoname}/sorted/{self.markerFile.split("/")[-1]}'), msg='Markerfile not copied')
        self.assertTrue(os.path.exists(f'{self.tomoname}/sorted/sorted_00.mrc'), msg='files not extracted')
        self.assertTrue(os.path.exists(self.metafile), msg='no metafile created')

    def test_03_MotionCorrection(self):
        pass

    def test_04_Alignment(self):
        """
        check that resulting alignment score is smaller than reference score, and alignmentResults are similar?
        """
        from pytom.basic.datatypes import DATATYPE_ALIGNMENT_RESULTS
        from pytom.gui.guiFunctions import loadstar
        import os

        cmd = f'pytom {self.pytomDir2}/reconstruction/generateAlignedTiltImages.py '

        cmd += f'--tiltSeriesName {self.tomoname}/sorted/sorted '
        cmd += f'--markerFile {self.tomoname}/sorted/markerfile.txt '
        cmd += f'--tiltSeriesFormat mrc '
        cmd += f'--firstIndex {self.firstIndexFull} '
        cmd += f'--lastIndex {self.lastIndexFull} '
        cmd += f'--referenceIndex {self.referenceTiltIndex} '
        cmd += f"--referenceMarkerIndex {self.referenceMarker} "
        cmd += f'--weightingType {self.weightingTypeAlignment} '
        cmd += f'--projectionTargets {self.tomoname}/alignment/marker____{float(self.startAngleFull):.1f},{float(self.endAngleFull):.1f}/GlobalAlignment/sorted '
        cmd += f'--lowpassFilter 0.9 '
        cmd += f'--expectedRotationAngle {self.expectedRotationAngle} '
        cmd += f'--numberProcesses 1 '

        print(cmd)
        os.system(cmd)

        ref = loadstar(self.alignment_sorted_results, dtype=DATATYPE_ALIGNMENT_RESULTS)
        ali = loadstar(f'{self.tomoname}/alignment/marker_0004_-60.0,56.0/GlobalAlignment/sorted/alignmentResults.txt',
                       dtype=DATATYPE_ALIGNMENT_RESULTS)

        for i in range(59):
            self.assertTrue(abs(ref['AlignmentTransX'][i] - ali['AlignmentTransX'][i]) < 0.001)
            self.assertTrue(abs(ref['AlignmentTransY'][i] - ali['AlignmentTransY'][i]) < 0.001)
            self.assertTrue(abs(ref['InPlaneRotation'][i] - ali['InPlaneRotation'][i]) < 0.001)
            self.assertTrue(abs(ref['Magnification'][i] - ali['Magnification'][i]) < 0.001)

        pass

    def test_05_Reconstruction(self):
        """
        check that resulting alignment score is smaller than reference score
        """
        import os
        if self.dont: return

        cmd = f"pytom {self.pytomDir}/bin/reconstructTomogram.py "
        cmd += f"--tiltSeriesName {self.tomoname}/sorted/sorted "
        cmd += f"--firstIndex {self.firstIndexFull} "
        cmd += f"--lastIndex {self.lastIndexFull} "
        cmd += f"--referenceIndex {self.referenceTiltIndex} "
        cmd += f"--markerFile {self.tomoname}/sorted/markerfile.txt "
        cmd += f"--referenceMarkerIndex {self.referenceMarker} "
        cmd += f"--expectedRotationAngle {self.expectedRotationAngle} "
        cmd += f"--projectionTargets {self.tomoname}/reconstruction/WBP/temp_files_unweighted/sorted_aligned "
        cmd += f"--projectionBinning 8 "
        cmd += f"--lowpassFilter 0.9  "
        cmd += f"--weightingType {self.weightingTypeRecon}  "
        cmd += f"--tomogramFile {self.tomoname}/reconstruction/WBP/tomogram_000_WBP.mrc "
        cmd += f"--tiltSeriesFormat mrc "
        cmd += f"--fileType mrc  "
        cmd += f"--tomogramSizeX 464 "
        cmd += f"--tomogramSizeY 464 "
        cmd += f"--tomogramSizeZ 464 "
        cmd += f"--reconstructionCenterX 0 "
        cmd += f"--reconstructionCenterY 0 "
        cmd += f"--reconstructionCenterZ 0 "

        print(cmd)

        os.system(cmd)

        self.assertTrue(os.path.exists(f'{self.tomoname}/reconstruction/WBP/tomogram_000_WBP.mrc'))

        pass

    def test_06_TemplateMatching_CPU(self):
        """
        check that resulting sum of correlation scores are larger than ref value. Check locations? Check if correct handedness has a higher score
        """
        import os
        if self.dont: return

        if 'gpu ' in device:
            return

        if not os.path.exists(f'{self.projectname}/04_Particle_Picking/Template_Matching/cross_correlation/tomogram_000_WBP'):
            os.mkdir(f'{self.projectname}/04_Particle_Picking/Template_Matching/cross_correlation/tomogram_000_WBP')

        for suffix in ('', '_Mirrored'):
            jobFile = f'''<JobDescription Destination="{self.projectname}/04_Particle_Picking/Template_Matching/cross_correlation/tomogram_000_WBP">
      <Volume Filename="{self.tomoname}/reconstruction/WBP/tomogram_000_WBP.mrc" Subregion=" 0,0,197,464,464,68 "/>
      <Reference Weighting="" File="{self.refDataDir}/human_ribo_21A{suffix}.em"/>
      <Mask Filename="{self.refDataDir}/human_ribo_mask_32_8_5{suffix}.mrc" Binning="1" isSphere="True"/>
      <WedgeInfo Angle1="30.0" Angle2="30.0" CutoffRadius="0.0" TiltAxis="custom">
        <Rotation Z1="0.0" Z2="0.0" X="0.0"/>
      </WedgeInfo>
      <Angles Type="FromEMFile" File="angles_12.85_7112.em"/>
      <Score Type="FLCFScore" Value="-100000000">
        <DistanceFunction Deviation="0.0" Mean="0.0" Filename=""/>
      </Score>
    </JobDescription>'''

            fname = f'{self.projectname}/04_Particle_Picking/Template_Matching/cross_correlation/tomogram_000_WBP/job{suffix}.xml'
            d = open(fname, 'w')
            d.write(jobFile)
            d.close()
            cmd = f"mpiexec --tag-output -n 16 pytom {self.pytomDir}/bin/localization.py -j {fname} -x 4 -y 4 -z 1"
            print(cmd)

            os.system(cmd)

        pass

    def test_06_TemplateMatching_GPU(self):
        """
        check that resulting sum of correlation scores are larger than ref value. Check locations? Check if correct handedness has a higher score
        """
        import os

        if 'cpu' in device: return

        if not os.path.exists(f'{self.projectname}/04_Particle_Picking/Template_Matching/cross_correlation/tomogram_000_WBP'):
            os.mkdir(f'{self.projectname}/04_Particle_Picking/Template_Matching/cross_correlation/tomogram_000_WBP')

        for suffix in ('', '_Mirrored'):
            jobFile = f'''<JobDescription Destination="{self.projectname}/04_Particle_Picking/Template_Matching/cross_correlation/tomogram_000_WBP">
      <Volume Filename="{self.tomoname}/reconstruction/WBP/tomogram_000_WBP.mrc" Subregion=" 0,0,197,464,464,68 "/>
      <Reference Weighting="" File="{self.refDataDir}/human_ribo_21A{suffix}.em"/>
      <Mask Filename="{self.refDataDir}/human_ribo_mask_32_8_5{suffix}.mrc" Binning="1" isSphere="True"/>
      <WedgeInfo Angle1="30.0" Angle2="30.0" CutoffRadius="0.0" TiltAxis="custom">
        <Rotation Z1="0.0" Z2="0.0" X="0.0"/>
      </WedgeInfo>
      <Angles Type="FromEMFile" File="angles_12.85_7112.em"/>
      <Score Type="FLCFScore" Value="-100000000">
        <DistanceFunction Deviation="0.0" Mean="0.0" Filename=""/>
      </Score>
    </JobDescription>'''

            fname = f'{self.projectname}/04_Particle_Picking/Template_Matching/cross_correlation/tomogram_000_WBP/job{suffix}.xml'
            d = open(fname, 'w')
            d.write(jobFile)
            d.close()
            cmd = f"mpiexec --tag-output -n 1 pytom {self.pytomDir}/bin/localization.py -j {fname} -x 1 -y 1 -z 1 -g 0"
            os.system(cmd)

        pass

    def test_07_TemplateMatchingCandidateExtractionGPU(self):
        """
        check that resulting scores and angle list is similar to gpu standalone
        """
        pass

    def test_08_CandidateExtraction(self):
        """
        check that resulting scores and angle list is similar to gpu standalone
        """

        for suffix in ('', '_Mirrored'):
            cmd = f'pytom {self.pytomDir}/bin/extractCandidates.py '
            cmd += f'--jobFile {self.tmfolder}/job{suffix}.xml '
            cmd += f'--result {self.tmfolder}/scores_human_ribo_21A{suffix}.em '
            cmd += f'--orientation {self.tmfolder}/angles_human_ribo_21A{suffix}.em '
            cmd += f'--particleList {self.projectname}/04_Particle_Picking/Picked_Particles/particleList_TM_tomogram_000_WBP{suffix}.xml '
            cmd += f'--particlePath Subtomograms/particleList_TM_tomogram_000_WBP '
            cmd += f'--size 8 '
            cmd += f'--numberCandidates 1500 '
            cmd += f'--minimalScoreValue 0.001 '
            print(cmd)
            os.system(cmd)



        pass

    def test_09_GenerateSubsetParticles(self):
        """
        check that resulting scores and angle list is similar to gpu standalone
        """
        from pytom.basic.structures import ParticleList
        import numpy as np



        particleListNormal = ParticleList()
        particleListNormal.fromXMLFile(f'{self.projectname}/04_Particle_Picking/Picked_Particles/particleList_TM_tomogram_000_WBP.xml')

        particleListMirror = ParticleList()
        particleListMirror.fromXMLFile(
            f'{self.projectname}/04_Particle_Picking/Picked_Particles/particleList_TM_tomogram_000_WBP_Mirrored.xml')


        scores = []
        for listParticles in (particleListNormal, particleListMirror):
            scores.append([])
            for p in listParticles:
                scores[-1].append(float(p.getScore().getValue()))

        scoresNormal = np.array(scores[0])
        scoresMirror = np.array(scores[1])

        self.assertTrue(scoresNormal.sum() > scoresMirror.sum(), 'Wrong handedness of reconstruction.')
        cutoff = 20+((scoresNormal[20:] > scoresMirror[20:])).astype(np.int32).sum()
        #print('cutoff: ', cutoff)
        self.assertTrue(cutoff > 1160, 'Wrong handedness of reconstruction.')
        self.assertTrue(scoresNormal[cutoff] > 0.13, "Poor correlation score")

        pl = particleListNormal[:cutoff]
        pl.toXMLFile(self.plFilename)


        pass

    def test_10_SubtomogramExtractionBinned(self):
        """
        check that resulting scores and angle list is similar to gpu standalone
        """

        cmd = f'cd {self.projectname}/05_Subtomogram_Analysis; reconstructWB.py '
        cmd += f'--particleList {self.plFilename} '
        cmd += f'--projectionDirectory {self.tomoname}/alignment/marker_0004_-60.0,56.0/GlobalAlignment/sorted '
        cmd += f'--coordinateBinning 8 '
        cmd += f'--size 100 '
        cmd += f'--applyWeighting {self.weightingTypeSubtomoRecon} '
        cmd += f'--projBinning 2 '
        cmd += f'--recOffset 0,0,0 '
        cmd += f'--metafile {self.tomoname}/sorted/mixedCTEM_tomo3.meta '
        cmd += f'--numProcesses 10 '

        os.system(cmd)

    def test_11_GLocal_CPU(self):
        """
        check that resulting resolution is similar to reference resolution
        """
        if 'gpu' in device: return

        outdir = os.path.join(self.glocaldir, 'alignment_000_cpu')

        if not os.path.exists(outdir): os.mkdir(outdir)

        cmd = f'cd {outdir}; mpiexec -n {self.numcores} pytom {self.pytomDir}/bin/GLocalJob.py '
        cmd += f'--particleList /data2/emuser01/test/stuart/tomo/04_Particle_Picking/Picked_Particles/particleList_TM_tomogram_000_WBP_good.xml '
        cmd += f'--mask /data2/emuser01/test/stuart/tomo/05_Subtomogram_Analysis/Alignment/GLocal/Glocal_mask.mrc '
        cmd += f'--numberIterations 8 '
        cmd += f'--pixelSize 5.24 '
        cmd += f'--particleDiameter 300 '
        cmd += f'--binning 1 '
        cmd += f'--destination {outdir} '
        cmd += f'--SphericalMask '
        cmd += f'--angleShells 3 '
        cmd += f'--angleIncrement 3.00 '
        cmd += f'--jobName {outdir}/glocal_input_params_Single_Iter.xml'
        print(cmd)
        os.system(cmd)

    def test_12_GLocal_GPU(self):
        """
        check that resulting resolution  is better than reference resolution and similar to CPU
        """

        if 'cpu' in device: return

        outdir = os.path.join(self.glocaldir, 'alignment_000_gpu')
        self.outdir = outdir

        if not os.path.exists(outdir): os.mkdir(outdir)



        cmd = f'cd {self.projectname}/05_Subtomogram_Analysis; mpiexec -n 5 pytom {self.pytomDir}/bin/GLocalJob.py '
        cmd += f'--particleList {self.plFilename} '
        cmd += f'--mask {self.refDataDir}/Glocal_mask.mrc '
        cmd += f'--numberIterations 4 '
        cmd += f'--pixelSize 5.24 '
        cmd += f'--particleDiameter 300 '
        cmd += f'--binning 1 '
        cmd += f'--destination {self.outdir} '
        cmd += f'--SphericalMask '
        cmd += f'--angleShells 3 '
        cmd += f'--angleIncrement 3.00 '
        cmd += f'--jobName {outdir}/glocal_input_params_Single_Iter.xml '
        cmd += f'--gpuID 0,1,2,3'

        os.system(cmd)

        pass

    def test_13_CompareGLocalVersionsResolution(self):
        """
        check that resulting scores and angle list is similar to gpu standalone
        """
        pass

    def test_14_UpdateParticleList(self):
        """
        check that resulting scores and angle list is similar to gpu standalone
        """
        from pytom.bin.updateParticleList import updatePL
        from pytom.basic.structures import ParticleList

        updatePL([self.glocaldir + '/alignment_000_gpu/3-ParticleList.xml'], [self.plFilenameReduced], suffix='_reduced', wedgeangles=[70,70], multiplyshift=1,)

        self.assertTrue(os.path.exists(self.plFilenameReduced), 'Updated particle list has not been created.')

        pl = ParticleList()
        pl.fromXMLFile(self.plFilenameReduced)

        # pl2 = ParticleList()
        #
        # for p in pl:
        #     if p.getScore().getValue() > 0.09:
        #         pl2.append(p)
        # print('\n\n\n\n', len(pl2), '\n\n\n')
        # pl2.toXMLFile(self.plFilenameReduced)

        self.assertTrue(pl[0].getWedge().getWedgeAngle() == 70, 'Wedge angles not updated')

    def test_15_CTFCorrection(self):
        """
        check that resulting scores of alignment of corrected are similar to ref values
        """

        import glob
        folder = f'{self.tomoname}/ctf'
        # Create stack from sorted
        self.ctfFileNameStack = f'{folder}/sorted.st'
        files = sorted(glob.glob(f'{self.tomoname}/sorted/sorted_??.mrc'))

        infile, outfile = open(os.path.join(folder, 'filein.txt'), 'w'), open(os.path.join(folder, 'fileout.txt'), 'w')
        cmd = 'cd {}; newstack -filei filein.txt -fileo fileout.txt '.format(folder)
        outfile.write('{}\n{}\n{}\n'.format(1, self.ctfFileNameStack, len(files)))
        infile.write('{}\n'.format(len(files)))
        for fname in files:
            infile.write('{}\n0\n'.format(fname))
        infile.close()
        outfile.close()
        os.system(cmd)

        # CTF correction using default files

        outfolder = f'{self.tomoname}/ctf/sorted_ctf'
        if not os.path.exists(outfolder): os.mkdir(outfolder)




        cmd = f'''cd {self.tomoname}/ctf;

ctfphaseflip -inp {self.ctfFileNameStack} -o ctfCorrected.st -an {self.ctfFileNameAngle} -defF {self.ctfFileNameDefocus} \
-defT 200 -iW 15 -pi 0.262 -cs 2.70 -gpu {self.singleGpuID} \
-am 0.08 -vo 300 -AxisAngle 180

mrcs2mrc.py -f ctfCorrected.st -t {self.tomoname}/ctf/sorted_ctf -p sorted_ctf -o {self.tomoname}/sorted '''


        os.system(cmd)

        # Align the ctf corrected files
        cmd = f'pytom {self.pytomDir}/reconstruction/generateAlignedTiltImages.py '

        cmd += f'--tiltSeriesName {self.tomoname}/ctf/sorted_ctf/sorted_ctf '
        cmd += f'--markerFile {self.tomoname}/sorted/markerfile.txt '
        cmd += f'--tiltSeriesFormat mrc '
        cmd += f'--firstIndex 20 '
        cmd += f'--lastIndex 40 '
        cmd += f'--referenceIndex 30 '
        cmd += f"--referenceMarkerIndex 4 "
        cmd += f'--weightingType 0 '
        cmd += f'--projectionTargets {self.tomoname}/alignment/marker____-20.0,20.0/GlobalAlignment/sorted_ctf '
        cmd += f'--lowpassFilter 0.9 '
        cmd += f'--expectedRotationAngle 0 '
        cmd += f'--numberProcesses 1 '
        os.system(cmd)
        pass

    def test_16_Subtomogram_Extraction(self, binning=2):
        """
        check for completion. Check correlation average is similar to reference value
        """

        cmd = f'cd {self.projectname}/05_Subtomogram_Analysis; reconstructWB.py '
        cmd += f'--particleList {self.plFilenameReduced} '
        cmd += f'--projectionDirectory {self.tomoname}/alignment/marker_0004_-20.0,20.0/GlobalAlignment/sorted_ctf '
        cmd += f'--coordinateBinning 8 '
        cmd += f'--size {200//binning} '
        cmd += f'--applyWeighting {self.weightingTypeSubtomoRecon} '
        cmd += f'--projBinning {binning} '
        cmd += f'--recOffset 0,0,0 '
        cmd += f'--metafile {self.tomoname}/sorted/mixedCTEM_tomo3.meta '
        cmd += f'--numProcesses 10 '
        print(cmd)
        os.system(cmd)

    def test_17_GLocal_Reduced_Binned_CPU(self):
        """
        check resulting resolution
        """

        if 'gpu' in device: return

        outdir = os.path.join(self.glocaldir, 'alignment_002_gpu')
        self.outdir_ali2 = outdir

        if not os.path.exists(outdir): os.mkdir(outdir)

        cmd = f'cd {self.projectname}/05_Subtomogram_Analysis; mpiexec -n {self.numcores} pytom {self.pytomDir}/bin/GLocalJob.py '
        cmd += f'--particleList {self.plFilenameReduced} '
        cmd += f'--mask {self.refDataDir}/Glocal_mask.mrc '
        cmd += f'--numberIterations 4 '
        cmd += f'--pixelSize 5.24 '
        cmd += f'--particleDiameter 300 '
        cmd += f'--binning 1 '
        cmd += f'--destination {self.outdir_ali2} '
        cmd += f'--SphericalMask '
        cmd += f'--angleShells 3 '
        cmd += f'--angleIncrement 3.00 '
        cmd += f'--jobName {self.outdir_ali2}/glocal_input_params_reduced.xml '
        cmd += f'--reference {os.path.join(self.glocaldir, "alignment_000_gpu")}/3-All.mrc '

        os.system(cmd)

        pass

    def test_18_GLocal_Reduced_Binned_GPU(self):
        """
        check that resulting resolution is below threshold and similar to CPU
        """

        if 'cpu' in device: return

        outdir = os.path.join(self.glocaldir, 'alignment_002_gpu')
        self.outdir_ali2 = outdir

        if not os.path.exists(outdir): os.mkdir(outdir)

        cmd = f'cd {self.projectname}/05_Subtomogram_Analysis; mpiexec -n 5 pytom {self.pytomDir}/bin/GLocalJob.py '
        cmd += f'--particleList {self.plFilenameReduced} '
        cmd += f'--mask {self.refDataDir}/Glocal_mask.mrc '
        cmd += f'--numberIterations 4 '
        cmd += f'--pixelSize 5.24 '
        cmd += f'--particleDiameter 300 '
        cmd += f'--binning 1 '
        cmd += f'--destination {self.outdir_ali2} '
        cmd += f'--SphericalMask '
        cmd += f'--angleShells 3 '
        cmd += f'--angleIncrement 3.00 '
        cmd += f'--jobName {self.outdir_ali2}/glocal_input_params_reduced.xml '
        cmd += f'--reference {os.path.join(self.glocaldir, "alignment_000_gpu")}/3-All.mrc '
        cmd += f'--gpuID 0,1,2,3 '

        os.system(cmd)

    def test_19_GLocal_Reduced_NonBinned_CPU(self):
            """
            check that resulting resolution is below threshold and similar to CPU
            """

            if 'gpu' in device: return

            from pytom.agnostic.transform import resize
            from pytom.agnostic.io import read, write

            cmd = f'cd {self.projectname}/05_Subtomogram_Analysis; reconstructWB.py '
            cmd += f'--particleList {self.plFilenameReduced} '
            cmd += f'--projectionDirectory {self.tomoname}/alignment/marker_0004_-20.0,20.0/GlobalAlignment/sorted_ctf '
            cmd += f'--coordinateBinning 8 '
            cmd += f'--size 200 '
            cmd += f'--applyWeighting {self.weightingTypeSubtomoRecon} '
            cmd += f'--projBinning 1 '
            cmd += f'--recOffset 0,0,0 '
            cmd += f'--metafile {self.tomoname}/sorted/mixedCTEM_tomo3.meta '
            cmd += f'--numProcesses 10 '
            print(cmd)
            os.system(cmd)

            outdir3 = os.path.join(self.glocaldir, 'alignment_003_gpu')
            if not os.path.exists(outdir3): os.mkdir(outdir3)

            outdir = os.path.join(self.glocaldir, 'alignment_002_gpu')
            self.outdir_ali2 = outdir

            v = read(f'{self.outdir_ali2}/3-All.mrc')

            r = resize(v, 2)
            reference = f'{outdir3}/resizedReference.mrc'
            write(reference, r)

            cmd = f'cd {self.projectname}/05_Subtomogram_Analysis; mpiexec -n {self.numcores} pytom {self.pytomDir}/bin/GLocalJob.py '
            cmd += f'--particleList {self.outdir_ali2}/3-ParticleList.xml '
            cmd += f'--mask {self.refDataDir}/Glocal_mask_200_75_5.mrc '
            cmd += f'--numberIterations 4 '
            cmd += f'--pixelSize 2.62 '
            cmd += f'--particleDiameter 300 '
            cmd += f'--binning 1 '
            cmd += f'--destination {outdir3} '
            cmd += f'--SphericalMask '
            cmd += f'--angleShells 3 '
            cmd += f'--angleIncrement 3.00 '
            cmd += f'--jobName {outdir3}/glocal_input_params_reduced.xml '
            cmd += f'--reference {reference} '


            os.system(cmd)

            from pytom.bin.gen_mask import gen_mask_fsc
            from pytom.agnostic.io import read, write

            model = sorted([os.path.join(outdir3, f) for f in os.listdir(outdir3) if
                            'average-FinalFiltered' in f and f.endswith('mrc')])[0]

            gen_mask_fsc(read(model), 4, f'{self.projectname}/05_Subtomogram_Analysis/Validation/maskFinalAverage.mrc',
                         1, 3)

            # alignedVolume = f'{outdir}/referrenceAligned3-All.mrc'
            # alignTwoVolumes(f'{outdir}/3-All.mrc', self.reference, outname=alignedVolume)

            cmd = f'''cd {self.projectname}/05_Subtomogram_Analysis/Validation

    fsc.py  '''
            cmd += f'--v1 {outdir3}/average-Final-Even.mrc  '
            cmd += f'--v2 {outdir3}/average-Final-Odd.mrc '
            cmd += f'--mask {self.projectname}/05_Subtomogram_Analysis/Validation/maskFinalAverage.mrc '
            cmd += f'--outputFolder {self.projectname}/05_Subtomogram_Analysis/Validation '
            cmd += f'--fsc 0.143 '
            cmd += f'--pixelsize 2.62 '
            cmd += f'--randomizePhases 0.000 '
            cmd += f'--combinedResolution '

            result = os.popen(cmd).read()

            resolution = float(result.split('Resolution determined for pixelsize :')[1].split()[-2])

            self.assertTrue(resolution < 10.,
                            'Final Resolution of the reconstruction is {resolution}. A resolution below 10. Angstrom is expected.')

    def test_20_GLocal_Reduced_NonBinned_GPU(self):
        """
        check that resulting resolution is below threshold and similar to CPU
        """

        if 'cpu' in device: return

        from pytom.agnostic.transform import resize
        from pytom.agnostic.io import read, write

        cmd = f'cd {self.projectname}/05_Subtomogram_Analysis; reconstructWB.py '
        cmd += f'--particleList {self.plFilenameReduced} '
        cmd += f'--projectionDirectory {self.tomoname}/alignment/marker_0004_-20.0,20.0/GlobalAlignment/sorted_ctf '
        cmd += f'--coordinateBinning 8 '
        cmd += f'--size 200 '
        cmd += f'--applyWeighting {self.weightingTypeSubtomoRecon} '
        cmd += f'--projBinning 1 '
        cmd += f'--recOffset 0,0,0 '
        cmd += f'--metafile {self.tomoname}/sorted/mixedCTEM_tomo3.meta '
        cmd += f'--numProcesses 10 '
        print(cmd)
        os.system(cmd)

        outdir3 = os.path.join(self.glocaldir, 'alignment_003_gpu')
        if not os.path.exists(outdir3): os.mkdir(outdir3)

        outdir = os.path.join(self.glocaldir, 'alignment_002_gpu')
        self.outdir_ali2 = outdir


        v =  read(f'{self.outdir_ali2}/3-All.mrc')

        r = resize(v,2)
        reference = f'{outdir3}/resizedReference.mrc'
        write(reference, r)


        cmd = f'cd {self.projectname}/05_Subtomogram_Analysis; mpiexec -n 5 pytom {self.pytomDir}/bin/GLocalJob.py '
        cmd += f'--particleList {self.outdir_ali2}/3-ParticleList.xml '
        cmd += f'--mask {self.refDataDir}/Glocal_mask_200_75_5.mrc '
        cmd += f'--numberIterations 4 '
        cmd += f'--pixelSize 2.62 '
        cmd += f'--particleDiameter 300 '
        cmd += f'--binning -1 '
        cmd += f'--destination {outdir3} '
        cmd += f'--SphericalMask '
        cmd += f'--angleShells 3 '
        cmd += f'--angleIncrement 3.00 '
        cmd += f'--jobName {outdir3}/glocal_input_params_reduced.xml '
        cmd += f'--reference {reference} '
        cmd += f'--gpuID 0,1,2,3 '

        os.system(cmd)


        from pytom.bin.gen_mask import  gen_mask_fsc
        from pytom.agnostic.io import read, write

        model = sorted([os.path.join(outdir3, f) for f in os.listdir(outdir3) if 'average-FinalFiltered' in f and f.endswith('mrc')])[0]

        gen_mask_fsc(read(model), 4, f'{self.projectname}/05_Subtomogram_Analysis/Validation/maskFinalAverage.mrc', 1, 3 )

        #alignedVolume = f'{outdir}/referrenceAligned3-All.mrc'
        #alignTwoVolumes(f'{outdir}/3-All.mrc', self.reference, outname=alignedVolume)

        cmd = f'''cd {self.projectname}/05_Subtomogram_Analysis/Validation

fsc.py  '''
        cmd += f'--v1 {outdir3}/average-Final-Even.mrc  '
        cmd += f'--v2 {outdir3}/average-Final-Odd.mrc '
        cmd += f'--mask {self.projectname}/05_Subtomogram_Analysis/Validation/maskFinalAverage.mrc '
        cmd += f'--outputFolder {self.projectname}/05_Subtomogram_Analysis/Validation '
        cmd += f'--fsc 0.143 '
        cmd += f'--pixelsize 2.62 '
        cmd += f'--randomizePhases 0.000 '
        cmd += f'--combinedResolution '
        cmd += f'--gpuID {self.singleGpuID}'

        print(cmd)
        result = os.popen(cmd).read()

        resolution = float(result.split('Resolution determined for pixelsize :')[1].split()[-2])
        print(f'Determined resolution: {resolution:.2f}')

        self.assertTrue(resolution < 10., f'Final Resolution of the reconstruction is {resolution}. A resolution below 10. Angstrom is expected.')

    def test_21_CCC_CPU(self):
        """
        check that resulting ccc is better than reference
        """
        pass

    def test_22_CCC_GPU(self):
        """
        check that resulting ccc is better than reference
        """
        pass

    def test_23_CPCA(self):
        """
        check that resulting class averages are similar to references
        """
        pass

    def test_24_AC3D(self):
        """
        check for completion
        """
        pass

    def test_25_AC3D_GPU(self):
        """
        check that results are similar to reference
        """
        pass

    def test_26_CleanUp(self):
        import os
        self.cleanUp()
        #self.assertTrue( not os.path.exists(self.projectname), msg="folder " + self.projectname + " has not been removed.")
        pass


if __name__ == '__main__':
    unittest.main()