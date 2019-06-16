rm -rf build
mkdir -p build 
cd build
cmake ../ -DCMAKE_INSTALL_PREFIX=$PREFIX -DCMAKE_BUILD_TYPE=Release
make 
make install

cd ..
pip install .

cp smc2 $PREFIX/bin


