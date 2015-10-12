#!/bin/bash
echo $LANG
echo $LC_ALL
if [ $TRAVIS_OS_NAME == linux ]; then
    sudo apt-get update -qq
    sudo apt-get install -qq git-core libcppunit-dev graphviz valgrind r-base-core

    if [ $CXX == g++ ]; then
        sudo apt-get install -qq g++-4.8
        export CXX="g++-4.8" CC="gcc-4.8"
    fi

fi

if [ $TRAVIS_OS_NAME == osx ]; then
    brew update && brew install llvm cppunit graphviz valgrind && brew link --force llvm
fi
