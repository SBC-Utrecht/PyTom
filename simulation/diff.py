import sys
from pytom.tompy.io import read


data1 = read(sys.argv[1]).squeeze()
data2 = read(sys.argv[2]).squeeze()

print(data1.shape, data2.shape)

import matplotlib
matplotlib.use('Qt5Agg')

diff = abs(data1-data2)



from pylab import *

print(allclose(data2, data1, rtol=1E-5))


fig,ax = subplots(1,3,figsize=(15,5))

ax[0].imshow(data1)
ax[1].imshow(data2)
ax[2].imshow(abs(data1-data2))
show()
