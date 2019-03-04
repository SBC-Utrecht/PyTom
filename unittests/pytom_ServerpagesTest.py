'''
Created on Dec 7, 2011

@author: hrabe
'''
import unittest

class pytom_ServerpagesTest(unittest.TestCase):
    
    
    def ReconstructionCall_Test(self):
        from pytom.frontend.serverpages import createReconstructionJob
        import os
        
        #reconstruction of tomogram
        requestString = 'prlDIR=./testData/reconstruction/test/projections/&tomo=tmp.em&x=2&y=2&z=2&sb=4&sa=1&jobFile=./test.sh' 
        result = createReconstructionJob.run(requestString)
        if result[-1] == '\n':
            result = result[:-1]
            
        assert result == '<File status="created" path="./test.sh" type="ReconstructionJob"/>'
        
        os.remove('./test.sh')
        
        #reconstruction of particles 
        requestString = 'prlDIR=./testData/reconstruction/test/projections/&plXML=./testData/xmlFiles/particleList.xml&x=100&y=100&z=100&sb=4&sa=1&xc=100&yc=100&zc=100&ts=4&jobFile=./test.sh'
        
        result = createReconstructionJob.run(requestString)
        if result[-1] == '\n':
            result = result[:-1]
            
        assert result == '<File status="created" path="./test.sh" type="ReconstructionJob"/>'
        
        os.remove('./test.sh')
        
    def LocalizationCall_Test(self):
        from pytom.frontend.serverpages import createLocalizationJob
        import os

        requestString = 'tomo=./testData/ribo.em&ref=./testData/ribo.em&mask=./testData/ribo.em&angle=angles_18_3040.em&dest=./results&low=0&high=10&smooth=0&w1=30&w2=30&x=2&y=2&z=2&jobFile=./test.xml'

        result = createLocalizationJob.run(requestString)
        
        if result[-1] == '\n':
            result = result[:-1]
            
        assert result == '<File status="created" path="./test.xml" type="LocalizationJob"/>'
        
        os.remove('./test.xml')
        os.remove('./test.sh')
        
        
    def AlignmentCall_Test(self):
        from pytom.frontend.serverpages import createAlignmentJob
        import os

        requestString = 'plXML=./testData/xmlFiles/particleList.xml&pixSize=4.7&partDia=250&ref=./testData/ribo.em&mask=./testData/ribo.em&sampling=GLOBAL&angFile=angles_35.76_320.em&lowestF=0&highestF=10&filtSm=0&adapt=ON&adResC=0.5&adResOf=0.1&angFac=0.5&score=flcf&pkPriRad=-1&pkSmooth=-1&iter=10&binning=1&dest=./res/&jobFile=./test.xml'
        
        result = createAlignmentJob.run(requestString)
        
        if result[-1] == '\n':
            result = result[:-1]
            
        assert result == '<File status="created" path="./test.xml" type="AlignmentJob"/>'
        
        os.remove('./test.xml')
        os.remove('./test.sh')

    def MCOEXMXCall_Test(self):
        from pytom.frontend.serverpages import createMCOEXMXJob
        import os
        
        requestString = 'plXML=./testData/xmlFiles/particleList.xml&pixSize=4.7&partDia=250&mask=./testData/ribo.em&lowestF=0&highestF=10&filtSm=0&score=flcf&iter=10&classes=2&binning=1&dest=./results/&conv=0&wa1=30&wa2=30&jobFile=./test.xml'
        result = createMCOEXMXJob.run(requestString)
        
        if result[-1] == '\n':
            result = result[:-1]
            
        assert result == '<File status="created" path="./test.xml" type="MCOEXMXJob"/>'
        
        os.remove('./test.xml')
        os.remove('./test.sh')
        
    def MCOACCall_Test(self):
        from pytom.frontend.serverpages import createMCOACJob
        import os
        
        requestString = 'plXML=./testData/xmlFiles/particleList.xml&pixSize=4.7&partDia=250&mask=./testData/ribo.emm&lowestF=0&highestF=10&filtSm=0&score=flcf&classes=2&binning=1&dest=./results/&conv=0&wa1=30&wa2=30&temp=sigma&stemp=1&astep=1&crit=metropolis&refin=5&jobFile=./test.xml'
        
        result = createMCOACJob.run(requestString)
        
        if result[-1] == '\n':
            result = result[:-1]
            
        assert result == '<File status="created" path="./test.xml" type="MCOACJob"/>'
        
        os.remove('./test.xml')
        os.remove('./test.sh')
        
    