"""
check scripts in bin/ directory
"""
from pytom.tools.files import checkDirExists,checkFileExists

import unittest

class pytom_BinTest(unittest.TestCase):

	def localization_Test(self):

		import os
		from pytom.tools.files import getPytomPath
		if not os.path.exists('./results'): os.mkdir('./results')
		command = getPytomPath() + '/bin/pytom ' + getPytomPath() + \
                    '/bin/localizationJob.py -v testData/ribo.em -r testData/ribo.em \
                    -m testData/ribo.em --wedge1 30 --wedge2 30 -a angles_18_3040.em \
                    -d ./results -b 10 --splitX 2 --splitY 2 --splitZ 2 -j ./test.xml > /dev/null'

		os.system(command)

		fileExists = checkFileExists('./test.xml')

		os.rmdir('./results')
		os.remove('./test.xml')
		os.remove('./test.sh')

		self.assertTrue(fileExists, 'localization_Test: files not written')

	def alignment_Test(self):
		import os
		from pytom.tools.files import getPytomPath
		if not os.path.exists('./results'): os.mkdir('./results')
		command = getPytomPath() + '/bin/pytom ' + getPytomPath() + '/bin/alignJob.py -p ./testData/xmlFiles/particleList.xml -r testData/ribo.em -m testData/ribo.em -l 0 -h 10 --angleShells 1 --angleIncrement 10 -d ./results -b 10 -n 10 --pixelSize 4.7 --particleDiameter 250 -j ./test.xml > /dev/null'
		
		os.system(command)

		fileExists = checkFileExists('./test.xml')

		os.rmdir('./results')
		os.remove('./test.xml')
		os.remove('./test.sh')

		self.assertTrue(fileExists, 'alignment_Test: files not written')
		
	def mcoEXMX_Test(self):

		import os
		from pytom.tools.files import getPytomPath
		if not os.path.exists('./results'): os.mkdir('./results')
		command = getPytomPath() + '/bin/pytom ' + getPytomPath() + \
                    '/bin/mcoEXMXJob.py -p ./testData/xmlFiles/particleList.xml \
                    -m testData/ribo.em -c 2 -t 0.1 --wedge1 30 --wedge2 30 -l 0 -h 10 \
                    -d ./results/ -n 10 -b 1 --pixelSize 4.7 --particleDiameter 250 \
                    -j test.xml > /dev/null'
		os.system(command)

		fileExists = checkFileExists('./test.xml')

		os.rmdir('./results')
		os.remove('./test.xml')
		os.remove('./test.sh')

		self.assertTrue(fileExists, 'mcoEXMX_Test: files not written')

	def mcoAC_Test(self):

		import os
		from pytom.tools.files import getPytomPath
		if not os.path.exists('./results'): os.mkdir('./results')
		command = getPytomPath() + '/bin/pytom ' + getPytomPath() + \
                    '/bin/mcoACJob.py -p ./testData/xmlFiles/particleList.xml \
                    -m testData/ribo.em -c 2 -t 0.1 --wedge1 30 --wedge2 30 -l 0 -h 10 \
                    -d ./results/ -n 10 -b 1 --pixelSize 4.7 --particleDiameter 250 \
                    --annealingStep 0.1 --startTemperature 5 -j test.xml > /dev/null'

		os.system(command)

		fileExists = checkFileExists('./test.xml')

		os.rmdir('./results')
		os.remove('./test.xml')
		os.remove('./test.sh')

		self.assertTrue(fileExists, 'mcoAC_Test: files not written')



