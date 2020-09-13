#!/bin/bash

if [ ! -d external ]; then
  mkdir external
fi

if [ ! -d external/fftw-3.3.8 ]; then
  pushd external
  wget http://www.fftw.org/fftw-3.3.8.tar.gz
  tar xzf fftw-3.3.8.tar.gz
  mkdir fftw_inst
  cd fftw-3.3.8
  ./configure --prefix $( readlink -f ../fftw_inst ) --enable-sse2 --enable-avx --enable-avx2 --enable-avx512 --enable-generic-simd128 --enable-generic-simd256 --enable-fma --disable-fortran
  make
  make install
  popd
fi

g++ -O3 -std=c++14 -o advection1d -Iexternal/fftw_inst/include -Lexternal/fftw_inst/lib *.cxx models/*.cxx solvers/*.cxx -lblas -lfftw3
#g++ -O0 -g -std=c++14 -o advection1d -Iexternal/fftw_inst/include -Lexternal/fftw_inst/lib *.cxx models/*.cxx solvers/*.cxx -lblas -lfftw3

exit 0
