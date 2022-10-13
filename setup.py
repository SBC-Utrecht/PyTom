from setuptools import setup, find_packages
from setuptools.command.install import install

import subprocess
import sys
import os
from pytom import __version__


def find_angle_lists(folder):
    return [e for e in os.listdir(folder) if e.endswith('.em') and e.startswith('angles_')]


def find_executables():
    folder = 'pytom/bin'
    a = [f'{folder}/{e}' for e in os.listdir(folder) if os.path.isfile(f'{folder}/{e}') and not '__' in e] + \
        ['pytom/bin/pytom', 'pytom/bin/ipytom', 'pytom/bin/pytomGUI']
    return a


class CustomInstall(install):
    def run(self):
        condadir = self.prefix
        version = f'{sys.version_info[0]}.{sys.version_info[1]}'
        commandInstall = f'python{version} compile.py --target all'
        if os.path.exists(condadir): commandInstall += f' --minicondaEnvDir {condadir}' 
        process = subprocess.Popen(commandInstall, shell=True, cwd="pytom/pytomc")
        process.wait()
        install.run(self)


setup(
    name='pytom',
    version=__version__,
    packages=find_packages(),
    package_dir={'pytom':'pytom'},
    package_data={'pytom/angles/angleLists': find_angle_lists('pytom/angles/angleLists')},
    include_package_data=True,
    author='`FridoF',
    author_email='gijsschot@gmail.com',
    url='https://github.com/FridoF/PyTomPrivate.git',
    install_requires=['lxml', 'PyFFTW', 'scipy', 'boost', 'numpy'],
    extras_require={
        'gpu': ['cupy'],
        'gui': ['PyQt5', 'pyqtgraph', 'mrcfile'],
        'all': ['cupy', 'PyQt5', 'pyqtgraph', 'mrcfile']},
    cmdclass={'install': CustomInstall},
    scripts=find_executables())

