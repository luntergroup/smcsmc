#!/bin/bash
git submodule update --init --recursive
cmake .
make

