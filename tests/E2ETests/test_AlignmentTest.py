"""
Alignment and reconstruction test between cpu, cpu parallel, and gpu. Ensures reconstruction is consistent across
implementations.

@author: McHaillet
"""
# gpu env needs to be initialized on top, otherwise environment is not set correctly
from pytom.gpu.initialize import initialize_gpu
initialize_gpu(0)

import unittest
import os
from pytom.reconstruction.reconstructionStructures import ProjectionList
from pytom.agnostic.io import write


class pytom_AlignmentTest(unittest.TestCase):
    def setUp(self):
        # set cpu/gpu (gpu by default 0)
        self.cpu_cores = 8

        # initialize projection list
        self.projections = ProjectionList()

        # TODO replace this by github test folder
        os.chdir('/data2/mchaillet')
        # alignment data is in local directory and should not be there (should be example data for unittests
        self.projections.loadDirectory(
            'processing/pytom_projects/tutorial/03_Tomographic_Reconstruction/tomogram_000/sorted',
            metafile='processing/pytom_projects/tutorial/03_Tomographic_Reconstruction/tomogram_000/'
                     'sorted/mixedCTEM_tomo1.meta',
            prefix='sorted_')
        self.projections.load_alignment(
            'processing/pytom_projects/tutorial/03_Tomographic_Reconstruction/tomogram_000/alignment/'
            'marker_0000_-60.0,54.0/GlobalAlignment/sorted/alignmentResults.txt')

    def test_reconstruction(self):
        from pytom.agnostic.correlation import xcc, nxcc
        from pytom_numpy import vol2npy

        # reconstruct
        recon_gpu = self.projections.reconstructVolumeGPU(weighting=-1, binning=8)
        # write('processing/pytom_projects/tutorial/alignment_test/tomo_000_recon_gpu_w-1_b8.mrc', recon_gpu)

        recon_cpu = self.projections.reconstructVolumeCPU(weighting=-1, binning=8, num_procs=self.cpu_cores)

        # convert to numpy and correlate
        rec_gpu_np = recon_gpu.get()
        rec_cpu_np = vol2npy(recon_cpu).copy()

        # calculate correlation
        ccc = nxcc(rec_gpu_np, rec_cpu_np)  # this should be better than 0.999
        self.assertGreater(ccc, 0.995, msg='correlation between cpu and gpu reconstruction not sufficient')

    def test_alignment(self):
        from pytom_numpy import vol2npy
        from pytom.agnostic.correlation import nxcc

        # run the alignment
        output_gpu = self.projections.to_projection_stack_gpu(weighting=-1, binning=4, show_progress_bar=True,
                                                          verbose=False)
        # write('processing/pytom_projects/tutorial/alignment_test/tomo_000_gpu.mrc', output_gpu[0])

        output_cpu_parallel = self.projections.to_projection_stack_parallel(weighting=-1, binning=4,
                                                                         show_progress_bar=True,
                                                           num_procs=10, verbose=False)
        output_cpu = self.projections.to_projection_stack(weighting=-1, binning=4, show_progress_bar=True,
                                                          verbose=False)

        # convert all to numpy
        stack_s, phi_s, theta_s, offset_s = (vol2npy(o).copy() for o in output_cpu)  # s = single
        stack_p, phi_p, theta_p, offset_p = (vol2npy(o).copy() for o in output_cpu_parallel)  # p = parallel
        stack_g, phi_g, theta_g, offset_g = (o.get() for o in output_gpu)  # g =  gpu
        offset_s = offset_s.squeeze().T

        # cpu and cpu parallel need to be identical
        self.assertTrue((stack_s != stack_p).sum() == 0, msg='alignment not consistent between cpu single core and '
                                                             'cpu parallel')

        # check if stacks are identical
        self.assertTrue((phi_s != phi_g).sum() == 0, msg='phi stack not identical between cpu and gpu')
        self.assertTrue((theta_s != theta_p).sum() == 0, msg='theta stack not identical between cpu and gpu')
        self.assertTrue((offset_s != offset_g).sum() == 0, msg='theta stack not identical between cpu and gpu')

        # check the nxcc between stack cpu and stack gpu
        ccc = nxcc(stack_s, stack_g)  # this should be better than 0.999
        self.assertGreater(ccc, 0.999, msg='correlation between cpu and gpu alignment not sufficient')

    def test_exact_filter(self):
        from pytom_numpy import vol2npy
        from pytom.agnostic.correlation import nxcc

        recon_gpu = self.projections.reconstructVolumeGPU(weighting=1, binning=8)
        write('processing/pytom_projects/tutorial/alignment_test/tomo_000_recon_gpu_w1_b8.mrc', recon_gpu)

        recon_gpu = self.projections.reconstructVolumeGPU(weighting=-1, binning=8)
        write('processing/pytom_projects/tutorial/alignment_test/tomo_000_recon_gpu_w-1_b8.mrc', recon_gpu)

        # run the alignment
        # output_gpu = self.projections.to_projection_stack_gpu(weighting=1, binning=4, show_progress_bar=True,
        #                                                       verbose=False)
        # write('processing/pytom_projects/tutorial/alignment_test/tomo_000_gpu.mrc', output_gpu[0])
        #
        # output_cpu_parallel = self.projections.to_projection_stack_parallel(weighting=1, binning=4,
        #                                                                     show_progress_bar=True,
        #                                                                     num_procs=10, verbose=False)
        #
        # # convert all to numpy
        # stack_p, phi_p, theta_p, offset_p = (vol2npy(o).copy() for o in output_cpu_parallel)  # p = parallel
        # stack_g, phi_g, theta_g, offset_g = (o.get() for o in output_gpu)  # g =  gpu
        #
        # # check the nxcc between stack cpu and stack gpu
        # ccc = nxcc(stack_p, stack_g)  # this should be better than 0.999
        # self.assertGreater(ccc, 0.999, msg='correlation between cpu and gpu alignment not sufficient')


if __name__ == '__main__':
    unittest.main()
