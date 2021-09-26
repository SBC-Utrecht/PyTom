# PytomGUI

PyTom is a toolbox developed for interpreting cryo electron tomography data. All steps from reconstruction, localization, alignment and classification are covered with standard and improved methods.

## Getting Started



### Prerequisites

PyTomGUI is designed to run on linux systems, but can also be installed on MacOSX. It requires the following software package to be installed:

```
# General packages 
- python (>= 3.7 )
- openmpi 
- fftw3
- gcc (version 5-7) 
- libxml2
- swig (>= 3.0.12)

# Python Packages
- numpy
- scipy
- boost
- lxml 
- mrcfile
- tqdm
- scikit-image
- matplotlib
- cupy (Optional for GPU acceleration)

# Optional Software Packages used by GUI
- PyQt5
- pyqtgraph
- motioncor2 ( >=1.2.1)
- imod (=4.10.25)

```

### Installing

Before you can install PyTom, you need to create an account on github. And to the github account you need to link a token for online identification. For more information click on the following link.

Furthermore, the software packages git needs to be installed. Git can be installed by sudo apt install git or yum install git. After git has been installed, run the following lines:

```
git clone git@github.com:FridoF/PyTomPrivate.git
cd PyTomPrivate
bash installMiniconda.sh
conda env create -f pytom_env.yml
```

Please remember the location where you decide to install conda (CONDA_INSTALL_DIR). 

```
conda activate pytom_env
python3.8 setup.py install --prefix [CONDA_INSTALL_DIR]/envs/pytom_env
```

### Installing on WSL2 (windows for linux subsystem)

PATH environment variable will contain windows directories in the /mnt folder of your linux subsystem. Pytom will have an issue with reading the PATH because of the white spaces ( arrrrg :( ). Make sure to update the PATH by removing the windows directories.

For WSL2 you need at least Windows 10 Pro with HyperV support.

to check if your machine can run of wsl, please type 'systeminfo' into the command prompt. The lines about HyperV should all return Yes.

To activate wsl, follow the following tutorial.
https://www.omgubuntu.co.uk/how-to-install-wsl2-on-windows-10


### Docker container

You can also use dockerfile to easily build pytom image on any platform.  
**Note:** Out of the box the image does not support GPU operations or pytomGUI.
```
git clone git@github.com:FridoF/PyTomPrivate.git && cd PyTomPrivate
docker build -t pytom .
```

Now that you have built an image, here are some examples of what you can do.
- Run ipytom: `docker run -it --rm pytom ipytom`  
- Run a script located on the host: `docker run -it --rm -v "/home/user/scripts:/hostfiles/" pytom pytom /hostfiles/some_script.py`

If you don't want to remove containers after the run, remove `--rm` flag.
Find info about how to user docker on [their docs](https://docs.docker.com/).


## Versioning

For the versions available, see the [tags on this repository]. 

## Authors

* **Gijs van der Schot** - *PyTomGUI* 
* **Thomas Hrabe**       - *PyTom* 
* **Yuxiang Chen**       - *PyTom*
* **Friedrich Forster**  - *PyTom* 

See also the list of [contributors] who participated in this project.

## License

Copyright (c) 2021

Utrecht University

http://pytom.sites.uu.nl

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
