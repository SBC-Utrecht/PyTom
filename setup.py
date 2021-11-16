from setuptools import setup, Extension, find_packages
from setuptools.command.install import install
import sys

import subprocess
#import setuptools.command.build_py

import os
from pytom import __version__


#print('Would you like to add a miniconda environment directory to your installation? Example: /home/user/miniconda3/envs/pytom_env')
#proc = subprocess.Popen("""conda env list | grep "*" | awk '{print $3}'""")#input('\nType Path: ')
#condadir = proc.communicate()[0]

condadir = os.popen("""conda env list | grep "*" | awk '{print $3}'""").read()[:-1]#input('\nType Path: ')

condadir = '' if not os.path.exists(condadir) else condadir

folder = 'pytom/angles/angleLists'
angleLists =  [e for e in os.listdir(folder) if e.endswith('.em') and e.startswith('angles_')]



def find_executables():
    folder = 'pytom/bin'
    a =  [f'{folder}/{e}' for e in os.listdir(folder) if os.path.isfile(f'{folder}/{e}') and not '__' in e] + ['pytom/bin/pytom', 'pytom/bin/ipytom', 'pytom/bin/pytomGUI']

    return a


class CustomInstall(install):
    def run(self):
        import sys
        # commandPullSubmodules = 'git submodule update --init --recursive'
        # process = subprocess.Popen(commandPullSubmodules, shell=True, cwd="./")
        # process.wait()
        #
        version = f'{sys.version_info[0]}.{sys.version_info[1]}'
        commandInstall = f'python{version} compile.py --target all'
        if os.path.exists(condadir): commandInstall += f' --minicondaEnvDir {condadir}' 
        process = subprocess.Popen(commandInstall, shell=True, cwd="pytom/pytomc")
        process.wait()

        install.run(self)



        # print('\nPlease add the following line to your .bashrc (LINUX) or to your .bash_profile (MAC): \n\n'
        #       f"\texport PATH='{os.getcwd()}/pytom/bin':$PATH\n\n")

        # t = input('ANSWER [y/n]: ')
        #
        # if t == 'y':
        #     for folder in (f"{os.environ['HOME']}/.bashrc", f"{os.environ['HOME']}/.bash_profile"):
        #         if os.path.exists(folder):
        #             a = open(folder, 'a')
        #             a.write(f'''\n#EXECUTABLES OF PYTOM ARE LOCATED IN FOLLOWING DIRECTORY\nexport PATH='{os.getcwd()}/pytom/bin':$PATH\n\n''')
        #             a.close()
        #             os.system(f'source {folder}')
        #             break

setup(
    name='pytom',
    version=__version__,
    packages=find_packages(),
    package_dir={'pytom':'pytom'},
    # package_data={'pytom':["alignment"]},
    package_data={'pytom/angles/angleLists': angleLists},
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

