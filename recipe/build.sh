$PYTHON -m pip install --ignore-installed --verbose .
rm -rf build
mkdir -p build 
cd build
cmake ../ -DCMAKE_INSTALL_PREFIX=$PREFIX -DCMAKE_BUILD_TYPE=Release -DBoost_NO_BOOST_CMAKE=ON 
make 
make install
