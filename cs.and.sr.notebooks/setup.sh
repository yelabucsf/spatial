#!/usr/bin/env bash
set -x
trap read debug

# this is a script to run on an amazon instance after git cloning
# downloads and installs all the necessary packages for running the notebooks

# download image
wget https://boygeniusreport.files.wordpress.com/2016/11/puppy-dog.jpg /myvol/data/

# install pip3
apt install python3-pip
pip3 install jupyterlab numpy scipy matplotlib cvxpy tqdm ipywidgets

# for tqdm_notebook to work in jupyter labs, need ipywidgets and the jupyter extension which requires nodejs and npm
apt install nodejs
apt install npm
jupyter labextension install @jupyter-widgets/jupyterlab-manager

# install libLBFGS
apt install libtool automake
wget -O /myvol/tools/liblbfgs-1.10.tar.gz https://github.com/downloads/chokkan/liblbfgs/liblbfgs-1.10.tar.gz
cd /myvol/tools
tar -xvzf liblbfgs-1.10.tar.gz
rm liblbfgs-1.10.tar.gz
cd liblbfgs-1.10/
./autogen.sh
./configure --enable-sse2
make
make install
cd ..
git clone https://rtaylor@bitbucket.org/rtaylor/pylbfgs.git
cd pylbfgs
ln -s /usr/local/lib/python3.6/dist-packages/numpy/core/include/numpy/ /usr/include/numpy # create a softlink to numpy_include folder, doesn't recognize it for some reason
python3 setup.py install

