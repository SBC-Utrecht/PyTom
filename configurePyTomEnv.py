import os
#sudo apt update
#sudo apt install git
#sudo apt install csh
#sudo apt install mlocate
#sudo apt install build-essential



conda = os.popen('which conda').read()[:-1]
print(conda)

if not conda:
  
   os.system(""" wget https://repo.anaconda.com/miniconda/Miniconda3-py38_4.10.3-Linux-x86_64.sh;  chmod +x Miniconda3-py38_4.10.3-Linux-x86_64.sh ; ./Miniconda3-py38_4.10.3-Linux-x86_64.sh""")

print('\nDo you want to create a conda environment from pytom_env.yml?')

c = ''
while c == '':
   c = input('\t[Y/n]: ')
   if not c in ('Y', 'y', 'n'):
      c = ''
      print('Please provide a valid option.')
if c in ('Y', 'y', 'yes' ):
   os.system('conda env create -f pytom_env.yml; conda init bash; conda activate pytom_env')

print('Do you want to install fftw3?')
install = ''
while not install:
    install = input('\t[Y/n]: ') 
    if not install in ('Y', 'n', 'y'):
       install = ''
       print('Please supply a valid option')
envdir = '/usr'
if conda:
   condapath = [line.split(':')[-1].replace('\t', '') for line in os.popen('conda info').read()[:-1].split('\n') if 'active env location : ' in line]
   try:
      suggestion = f"{condapath[0].strip(' ')}/envs/pytom_env"
      envdir = envdir if not os.path.exists(suggestion) else suggestion
   except Exception as e:
      print(e)
      pass
   
if install in ('Y', 'y'):
   installdir= ''
   print('\nWhere do you want to install fftw3? Press enter if selection is ok.')
   
   while installdir == '':
      installdir = input(f'\t{envdir}: ')
      if not installdir:
         installdir = envdir

      if not os.path.exists(installdir):
         installdir = ''

   if 'pytom_env' in installdir:
      if not os.path.exists(f'{installdir}/bin/x86_64-conda_cos6-linux-gnu-cc'): os.system(f'ln -s {installdir}/bin/cc {installdir}/bin/x86_64-conda_cos6-linux-gnu-cc')

   os.system(f' wget http://www.fftw.org/fftw-3.3.9.tar.gz; tar -zxvf fftw-3.3.9.tar.gz; cd fftw-3.3.9; ./configure --prefix={installdir} --enable-mpi --enable-shared=yes --libdir={installdir}; make install; make clean; ./configure --prefix={installdir} --enable-mpi --enable-shared=yes --enable-float; make install')
