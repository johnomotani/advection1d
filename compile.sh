#!/bin/bash

g++ -O3 -std=c++14 -o advection1d *.cxx models/*.cxx solvers/*.cxx -lblas
#g++ -O0 -g -std=c++14 -o advection1d *.cxx models/*.cxx solvers/*.cxx -lblas

exit 0
