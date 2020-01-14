#!/usr/bin/env pytom

from pytom.tompy.io import read
import sys
import numpy as np

a = read(sys.argv[1])
b = read(sys.argv[2])

d = np.abs(a-b)

print(d.mean(), d.min(), d.max())
