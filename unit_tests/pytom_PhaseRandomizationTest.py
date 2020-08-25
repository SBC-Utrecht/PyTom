'''
Author: GvdS
'''

import unittest

class pytom_PhaseRandomizationTest(unittest.TestCase):
    
    def test_phase_random_even(self):
        """
        test phase randomization
        """
        import numpy as np
        from pytom.tompy.correlation import generate_random_phases_3d

        a = np.zeros((128,128,128))
        a[10:-10,10:-10,10:-10] = 1
        ft = np.fft.rfftn(a)

        amplitude = abs(ft)
        phase = generate_random_phases_3d(ft.shape, reduced_complex=True)
        
        self.assertTrue( np.fft.irfftn((amplitude * np.exp(1j * phase))).imag.std() < 1E-8 )

    def test_phase_random_odd(self):
        """
        test phase randomization
        """
        import numpy as np
        from pytom.tompy.correlation import generate_random_phases_3d

        a = np.zeros((129,129,129))
        a[10:-10,10:-10,10:-10] = 1
        ft = np.fft.rfftn(a)

        amplitude = abs(ft)
        phase = generate_random_phases_3d(ft.shape, reduced_complex=True)
        
        self.assertTrue( np.fft.irfftn((amplitude * np.exp(1j * phase))).imag.std() < 1E-8 )

    def runTest(self):
        self.test_phase_random_odd()
        self.test_phase_random_even()
        
if __name__ == '__main__':
    unittest.main()
