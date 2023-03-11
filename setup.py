from setuptools import setup, find_namespace_packages
from setuptools.command.install import install
from setuptools.command.develop import develop

import subprocess
import sys
import pathlib
from pytom import __version__


def find_executables():
    bin_folder = pathlib.Path('pytom/bin')
    executables = [str(path) for path in bin_folder.iterdir() if path.is_file() and '__' not in path.name]
    return executables + [str(bin_folder.joinpath(f)) for f in ['pytom', 'ipytom', 'pytomGUI']]


def compile_pytom(miniconda_dir):
    install_command = f'{sys.executable} compile.py --target all'
    if miniconda_dir.exists():
        install_command += f' --minicondaEnvDir {miniconda_dir}'
    process = subprocess.Popen(install_command, shell=True, cwd="pytom/pytomc")
    process.wait()


class PyTomInstaller(install):
    def run(self):
        compile_pytom(pathlib.Path(self.prefix))
        install.run(self)


class PyTomDeveloper(develop):
    def run(self):
        compile_pytom(pathlib.Path(self.prefix))
        develop.run(self)

        """Hacky solution, but we cannot rely on the script linking of pip in develop-mode. Pip will put scripts in the 
        miniconda bin that have a different shebang than the original ones. However, we need the pytom shebang for GPU 
        scripts to execute. Only then is the GPU environment variable properly set. Updating gpu code to the cupy 
        advised agnostic setup (as put in PyTomPrivate issue #117) should also solve this."""
        miniconda_bin = pathlib.Path(self.prefix)
        scripts = find_executables()
        for script in scripts:
            script_path = pathlib.Path(script)
            symlink_path = miniconda_bin.joinpath('bin', script_path.name)
            symlink_path.unlink()
            symlink_path.symlink_to(script_path.absolute())


setup(
    name='pytom',
    version=__version__,
    packages=find_namespace_packages(include=['pytom*']),
    package_dir={'pytom': 'pytom'},
    package_data={
        'pytom.angles.angleLists': ['*.em'],
        'pytom.simulation.detectors': ['*.csv'],
        'pytom.simulation.membrane_models': ['*.pdb']
    },
    data_files=[("pytom_data", ["./LICENSE.txt"])],  # This is a relative dir to sys.prefix
    include_package_data=True,
    author='`FridoF',
    author_email='strubi.pytom@uu.nl',
    url='https://github.com/SBC-Utrecht/PyTom.git',
    install_requires=['lxml', 'scipy', 'boost', 'numpy'],
    extras_require={
        'gpu': ['cupy'],
        'gui': ['PyQt5', 'pyqtgraph', 'mrcfile'],
        'all': ['cupy', 'PyQt5', 'pyqtgraph', 'mrcfile']},
    cmdclass={'install': PyTomInstaller,
              'develop': PyTomDeveloper},
    scripts=find_executables())

