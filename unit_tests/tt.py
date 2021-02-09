from pytom.basic.correlation import nxcc
from pytom.basic.files import read


print(nxcc(read('GPU/vol_1_3.mrc'), read('GPU/vol_3_1.mrc'), read('CPU/mask.em')))
print(nxcc(read('CPU/vol_1_3.em'),  read('CPU/vol_3_1.em'),  read('CPU/mask.em')))


from pytom.tompy.correlation import nxcc
from pytom.tompy.io import read as readN


print(nxcc(readN('GPU/vol_1_3.mrc'), readN('GPU/vol_3_1.mrc'), readN('CPU/mask.em')))
print(nxcc(readN('CPU/vol_1_3.em'),  readN('CPU/vol_3_1.em'),  readN('CPU/mask.em')))
