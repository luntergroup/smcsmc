export CONDA_BUILD_SYSROOT=/opt/MacOSX10.9.sdk
export CFLAGS="${CFLAGS} -i sysroot ${CONDA_BUILD_SYSROOT}"
$PYTHON setup.py install
rm -rf build
mkdir -p build 
cd build
cmake ../ -DCMAKE_INSTALL_PREFIX=$PREFIX -DCMAKE_BUILD_TYPE=Release -DCMAKE_OSX_SYSROOT=${CONDA_BUILD_SYSROOT} -DCMAKE_OSX_DEPLOYMENT_TARGET=10.9
make 
make install
