'''
Created on Jul 30, 2012

@author: hrabe
'''

import unittest

class pytom_IOTest(unittest.TestCase):
    
    def test_read(self):
        """
        test general read function
        """
        from pytom.basic.files import read
        
        v = read('./testData/ribo.em')
        #v = read('testData/emd_1480.map.em_bin_4.em')
        self.assertTrue( v.sizeX() == 100 and v.sizeY() == 100 and v.sizeZ() == 100)
        
        v = read('./testData/ribo.em',binning=[4,4,4])
        assert v.sizeX() == 25 and v.sizeY() == 25 and v.sizeZ() == 25
    
    def test_readem(self):
        """
        test readem
        """
        from pytom.basic.files import read_em
        int2_file = 'testData/int2.em'
        (image, header) = read_em(int2_file)
        inibytes = header.get_1st4bytes()
        self.assertTrue( inibytes[3] == 5, 
            'Datatype number is wrong! should have been changed to 5 from 2 in initial file')
        dims = header.get_dim()

        #now binning
        (image, header) = read_em(int2_file, binning=[4, 4, 1])
        bdims = header.get_dim()
        self.assertTrue( bdims[0]*4 == dims[0], 
            'Dimensions of binned file are wrong!')
        
    def runTest(self):
        self.test_read()
        self.test_readem()
        
if __name__ == '__main__':
    unittest.main()
        
        
