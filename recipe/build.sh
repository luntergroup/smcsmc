$PYTHON setup.py install
rm -rf build
mkdir -p build 
cd build
cmake ../ -DCMAKE_INSTALL_PREFIX=$PREFIX -DCMAKE_BUILD_TYPE=Release
make 
make install
