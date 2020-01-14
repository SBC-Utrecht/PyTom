# PytomGUI

PyTom is a toolbox developed for interpreting cryo electron tomography data. All steps from reconstruction, localization, alignment and classification are covered with standard and improved methods.

## Getting Started



### Prerequisites

PyTomGUI is designed to run on linux systems, but can also be installed on MacOSX.It requires the following software package to be installed:

```
- python (>= 3.7 )
- openmpi 
- fftw3
- gcc (version 5-7) 
- numpy
- boost
- lxml, libxml2, libxstl. 
- swig (>= 3.0.12)
- PyQt5
- pyqtgraph
- motioncor2 ( >=1.2.1)
- imod (>=4.10.25)
- mrcfile
```

### Installing

To install PyTomGUI please clone the erice version by executing the following command 

```
git clone --recursive https://github.com/FridoF/PyTomPrivate.git pytom
```

After a succesful clone enter the new directory and go to pytomc
```
cd PyTomPrivate/pytomc
```

Here you will find an installation script name compile.py. Please run this script with python3. 
```
python3.7 compile.py --pythonVersion 3.7 --target all 
```

If certain dependencies are not found, chances are there you either have not installed them, or that the installation paths to the dependencies are not found automatically. 

It is possible to add installation paths manually by adding them after a respective flag (include dirs after --includeDir, lib dirs after --libDir and exe dirs after --exeDir). More than one path can be added per flag, where each path should be separated by a space. Note that the order of the paths can make a difference.

```
python3.7 compile.py --pythonVersion 3.7 --target all --exeDir [path_to_your_exe_dir1] [path_to_your_exe_dir2] --libDir [path_to_your_lib_dir1] [path_to_your_lib_dir2] --includeDir [path_to_your_exclude_dir1] [path_to_your_exclude_dir2] 
```

To locate missing dependencies try:

```
locate Python.h
```

This results in something similar to:
```
/usr/local/Cellar/python/3.7.2_1/Frameworks/Python.framework/Versions/3.7/include/python3.7m/Python.h
/usr/local/Cellar/python@2/2.7.15/Frameworks/Python.framework/Versions/2.7/include/python2.7/Python.h
```

Now update the includeDir flag to:

```
python3.7 compile.py --pythonVersion 3.7 --target all --includeDir /usr/local/Cellar/python/3.7.2_1/Frameworks/Python.framework/Versions/3.7/include/python3.7m/
```

## Versioning

For the versions available, see the [tags on this repository](https://gschot@bitbucket.org/gschot/pytom/tags). 

## Authors

* **Gijs van der Schot** - *PyTomGUI* 
* **Thomas Hrabe**       - *PyTom* 
* **Yuxiang Chen**       - *PyTom*
* **Friedrich Forster**  - *PyTom* 

See also the list of [contributors](https://gschot@bitbucket.org/gschot/pytom/contributors) who participated in this project.

## License

Copyright (c) 2019

Utrecht University

http://www.pytom.org

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

The complete license can be obtained from 
http://www.gnu.org/licenses/gpl-2.0.html.
