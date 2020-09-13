#!/bin/bash

g++ -O3 -std=c++14 -o advection1d *.cxx models/*.cxx solvers/*.cxx -lblas -lfftw3
#g++ -O0 -g -std=c++14 -o advection1d *.cxx models/*.cxx solvers/*.cxx -lblas -lfftw3

exit 0
