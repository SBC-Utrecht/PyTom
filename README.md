# PytomGUI

PyTom is a toolbox developed for interpreting cryo electron tomography data. All steps from reconstruction, localization, alignment and classification are covered with standard and improved methods.

## Getting Started

### Prerequisites

PyTomGUI is designed to run on linux systems. All required packages are managed via conda. For further info see the wiki [installation page](https://github.com/FridoF/PyTom/wiki/Installation).

If you want to make use of PyTomGUI interfacing with motioncor2 and IMOD's ctf correction, you will need to install these:

```
- motioncor2 ( >=1.2.1)
- imod (=4.10.25)
```

### Installing

Git needs to be installed to clone the repository. Git can be installed by sudo apt install git or yum install git. After git has been installed, clone the pytom repository and enter it:

```
mkdir pytom
git clone https://github.com/FridoF/PyTom.git pytom
cd pytom
```

Only in case you do not have the miniconda environment manager installed, run the following (otherwise continue with the next step):

```
bash installMiniconda.sh
```

Please remember the location where miniconda is installed [CONDA_INSTALL_DIR], because we need it later on. By default it will be installed in the homo directory ~/miniconda3.

Now we are ready to create the conda environment for pytom (solving all the dependencies might take a moment):

```
conda env create -f environments/pytom_py3.8_cu10.1.yaml --name pytom_env
```

Activate the environment and run the pytom installation scripts to compile the backend (will take ~5 minutes):

```
conda activate pytom_env
python3.8 setup.py install --prefix [CONDA_INSTALL_DIR]/envs/pytom_env
```

Now everything should be setup to start pytom. You can start the GUI by running the following:

```
pytomGUI
```

Whenever you want to open pytomGUI in a new shell, first run:

```
conda activate pytom_env
```

For more information on running pytom, see the [tutorial and wiki on github](https://github.com/FridoF/PyTom/wiki)

### Installing on WSL2 (windows for linux subsystem)

For WSL2 you need at least Windows 10 Pro with HyperV support. To check if your machine can run of wsl, please type 'systeminfo' into the command prompt. The lines about HyperV should all return Yes.

To activate wsl, see the following tutorial: https://www.omgubuntu.co.uk/how-to-install-wsl2-on-windows-10

You also need a custom setup of the CUDA toolkit. A GPU driver needs to be installed on windows, and the CUDA toolkit isntallation in WSL2 requires installation of special packages, see the nvidia docs for more info: https://docs.nvidia.com/cuda/wsl-user-guide/index.html

For installation you should be able to follow the instructions as above, however you need to update the .yaml environment file. From the dependencies remove: cudatoolkit, cudnn, cupy. Then add the following to the bottom row (tqdm + joblib only placed for illustration) of the file, where [VERSION] should be replaced with CUDA tookit version, e.g. 116 if you have 11.6 of the CUDA toolkit. 

```
...
  - tqdm=4.62.1
  - joblib=1.0.1
  - pip:
    - cupy-cuda[VERSION]
```

Afterwards you should be able to follow the installation steps as above.

### Docker container

You can also use dockerfile to easily build pytom image on any platform.  
**Note:** Out of the box the image does not support GPU operations or pytomGUI.
```
git clone https://github.com/FridoF/PyTom.git && cd PyTom
docker build -t pytom .
```

Now that you have built an image, here are some examples of what you can do.
- Run ipytom: `docker run -it --rm pytom ipytom`  
- Run a script located on the host: `docker run -it --rm -v "/home/user/scripts:/hostfiles/" pytom pytom /hostfiles/some_script.py`

If you don't want to remove containers after the run, remove `--rm` flag.
Find info about how to user docker on [their docs](https://docs.docker.com/).


### Unit tests

After installation you can run pytom unit tests that check the functionality of the program.

```
cd tests
pytom -m unittest discover
```

## Versioning

For the versions available, see the [tags on this repository]. 

## Authors

* **Marten Chaillet**    - *PyTomGUI* *GPU* *MicrographModeller*
* **Gijs van der Schot** - *PyTomGUI* *GPU*
* **Ilja Gubins**        - *GPU* *MicrographModeller*
* **Mihajlo Vanevic**    - *PyTomGUI*
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
