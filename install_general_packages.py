import sys
import os

if sys.argv[1] == '--ubuntu':
    cmd = ''' 
sudo apt install -y python3-pip git locate fftw3 fftw3-dev libxml2 python3-libxml2 csh libboost-dev 

pip install PyQt5==5.12.1 sip pyqtgraph matplotlib scikit-image

cd 
mkdir Programs
cd Programs
wget https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.0.tar.gz .

tar -xvzf openmpi-4.1.0.tar.gz

cd openmpi-4.1.0

./configure --prefix=’/usr/’ 
sudo make all install



''' 
