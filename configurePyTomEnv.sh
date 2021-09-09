
#sudo apt update
#sudo apt install git
#sudo apt install csh
#sudo apt install mlocate
#sudo apt install build-essential

wget https://repo.anaconda.com/miniconda/Miniconda3-py38_4.10.3-Linux-x86_64.sh
chmod +x Miniconda3-py38_4.10.3-Linux-x86_64.sh
./Miniconda3-py38_4.10.3-Linux-x86_64.sh

conda env create -f pytom_env.yml

conda activate pytom_env

#wget http://www.fftw.org/fftw-3.3.9.tar.gz
#tar -zxvf fftw-3.3.9.tar.gz
#cd fftw-3.3.9

#./configure --prefix=/home/gijs/miniconda3/pytom_env --enable-mpi --enable-shared=yes
#make install
