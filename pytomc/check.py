#!/usr/bin/env python

import os
import platform

if platform.system() == "Linux":
    os.system("./check.csh 'ldd -r -d swigModules/*so'")
    os.system("./check.csh 'readelf -d swigModules/*so'")
elif platform.system() == 'Darwin':
    os.system("./check.csh 'otool -L swigModules/*so'")
    


