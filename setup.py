from setuptools import setup, Extension, find_packages
from setuptools.command.install import install
import subprocess
import os


def find_executables():
    folder = 'pytom/bin'
    a =  [f'{folder}/{e}' for e in os.listdir(folder) if os.path.isfile(f'{folder}/{e}') and not '__' in e]
    return a

class CustomInstall(install):
    def run(self):
        commandPullSubmodules = 'git submodule update --recursive --remote; git submodule update --recursive'
        process = subprocess.Popen(commandPullSubmodules, shell=True, cwd="pytom")
        process.wait()

        commandInstall = 'python3.7 compile.py --target all --pythonVersion 3.7 > logfile.installation.txt'
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
    version='0.994',
    packages=find_packages(),
    package_dir={'pytom':'pytom'},
    # package_data={'pytom':["alignment"]},
    author='`FridoF',
    author_email='gijsschot@gmail.com',
    url='https://github.com/FridoF/PyTomPrivate.git',
    install_requires=['lxml', 'PyFFTW', 'scipy', 'boost', 'numpy'],
    extras_require={
        'gpu': ['cupy'],
        'gui': ['PyQt5', 'pyqtgraph', 'mrcfile'],
        'all': ['cupy', 'PyQt5', 'pyqtgraph', 'mrcfile']},
    cmdclass={'install': CustomInstall},
    include_package_data=True,
    scripts=find_executables(),
    test_suite='nose.collector',
    tests_require=['nose'])

