#!/bin/bash
echo $LANG
echo $LC_ALL
if [ $TRAVIS_OS_NAME = linux ]; then
    echo "Linux"
    sudo apt-get update -qq
    sudo apt-get install -qq git-core libcppunit-dev graphviz valgrind r-base-core
fi

if [ $CXX = "g++" ]; then
    echo "change g++"
    #sudo apt-get install -qq g++-4.8
    export CC="gcc-4.8";
    export CXX="g++-4.8";
    export LINK="gcc-4.8";
    export LINKXX="g++-4.8";
fi

if [ $TRAVIS_OS_NAME = osx ]; then
    echo "OSX"
    brew update && brew install llvm cppunit graphviz valgrind && brew link --force llvm
fi
