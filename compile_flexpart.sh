#!/bin/bash

module purge
module load gcc/7.3.0 intelmpi/18.0.2 mkl/18.0.2 hdf5-par/1.8.20 netcdf4/4.6.1

version=3.3.2
opt=mpi
#opt=omp
opt=serial
#wget https://www.flexpart.eu/downloads/58 -O flexpart-wrf-${version}.tar.gz
#tar xvf flexpart_wrf_${version}.tar.gz
#mv Src_flexwrf_v${version} Src_flexwrf_v${version}-${opt}
#cd Src_flexwrf_v${version}-${opt}
make -f makefile.mom clean
make -f makefile.mom $opt NETCDF=$(nc-config --prefix)