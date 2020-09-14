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
  # configuration for profiling
  #./configure --prefix $( readlink -f ../fftw_inst ) --enable-sse2 --enable-avx --enable-avx2 --enable-avx512 --enable-generic-simd128 --enable-generic-simd256 --enable-fma --disable-fortran CFLAGS=-pg
  make
  make install
  popd
fi

if [ ! -d external/BLAS-3.8.0 ]; then
  pushd external
  wget http://www.netlib.org/blas/blas-3.8.0.tgz
  tar xzf blas-3.8.0.tgz
  pushd BLAS-3.8.0
  make
  popd
  wget http://www.netlib.org/blas/blast-forum/cblas.tgz
  tar xzf cblas.tgz
  cd CBLAS
  cp Makefile.LINUX Makefile.in
  sed -i "s/libblas.a/..\/..\/BLAS-3.8.0\/blas_LINUX.a/" Makefile.in
  make
  popd
fi

# using system-installed libraries
#g++ -O3 -std=c++14 -o advection1d *.cxx models/*.cxx solvers/*.cxx -lblas -lfftw3

g++ -O3 -std=c++14 -o advection1d -Iexternal/fftw_inst/include -Lexternal/fftw_inst/lib -Iexternal/CBLAS/include *.cxx models/*.cxx solvers/*.cxx external/CBLAS/lib/cblas_LINUX.a external/BLAS-3.8.0/blas_LINUX.a -lgfortran -lfftw3

# for debugging
#g++ -O0 -g -std=c++14 -o advection1d -Iexternal/fftw_inst/include -Lexternal/fftw_inst/lib -Iexternal/CBLAS/include -Lexternal/CBLAS/lib -Lexternal/BLAS-3.8.0 *.cxx models/*.cxx solvers/*.cxx -lcblas_LINUX -lblas_LINUX -lgfortran -lfftw3

# for profiling
#g++ -pg -O3 -std=c++14 -o advection1d -Iexternal/fftw_inst/include -Lexternal/fftw_inst/lib -Iexternal/CBLAS/include *.cxx models/*.cxx solvers/*.cxx external/CBLAS/lib/cblas_LINUX.a external/BLAS-3.8.0/blas_LINUX.a -lgfortran -lfftw3

exit 0
