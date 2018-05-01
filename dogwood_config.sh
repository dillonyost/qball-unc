#!/bin/bash
# configure qball for ALCF #

LIBHOME=/nas/longleaf/apps/intel/17.2/intel/compilers_and_libraries_2017.2.174/linux/mkl/lib/intel64
MKLDIR=/nas/longleaf/apps/intel/17.2/intel/compilers_and_libraries_2017.2.174/linux/mkl
MPIDIR=/usr/mpi/intel/mvapich2-2.3a/
BLASDIR=$LIBHOME
LAPACKDIR=$LIBHOME
SCALAPACK_DIR=$LIBHOME
SCALAPACKLIB=$SCALAPACK_DIR/libmkl_scalapack_lp64.a
LAPACKLIB=$LAPACKDIR/libmkl_lapack95_lp64.a
#ESSLDIR=/soft/libraries/essl/current
#HPMLIBS="-L/soft/perftools/hpctw/lib -lmpihpm_smp -L/bgsys/drivers/ppcfloor/bgpm/lib -lbgpm"
#XERCESCDIR=/nas/longleaf/home/dyost/software/qbox/xerces-c-3.1.4/include
#XERCESCLIBDIR=/nas/longleaf/home/dyost/software/qbox/xerces-c-3.1.4/lib
#XERCESCDIR=/nas/longleaf/home/dyost/software/qbox/xerces-c-src_2_8_0
#XERCESCLIBDIR=/nas/longleaf/home/dyost/software/qbox/xerces-c-src_2_8_0/lib
#XERCESLIB=$XERCESCLIBDIR/libxerces-c.so

#BGQ_SDK_PATH=/bgsys/drivers/ppcfloor
export CXX=mpicc
export CC=mpicxx
export LIBS_BLAS="-L$LIBHOME/libmkl_blas95_lp64.a "
export MKLROOT="/nas/longleaf/apps/intel/17.2/intel/compilers_and_libraries_2017.2.174/linux/mkl/"

DFLAGS="-DUSE_MPI -DPRINTALL -DUSE_CSTDIO_LFS -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DUSE_XERCES -DXERCESC_3"

INCLUDE=" -I$MKLDIR/include -I$MPIDIR/include"

LIBPATH=" -L$LAPACKDIR -L$BLASDIR -L$MKLDIR/lib/intel64 -L$MPIDIR/lib64 "

#export LIBS="$SCALAPACKLIB -lmkl_intel_ilp64 -lmkl_lapack_ilp64 -mkl_sequential -lmkl_core -lirc -lifcore -lsvml -lxerces-c  -llapack -lpthread "
export PLIBS=" -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64"
export LIBS="-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lirc -lifcore -lsvml -lpthread $PLIBS $LAPACKLIB -lxerces-c" #-lxlf90_r -lxlfmath -lxlsmp"
export LDFLAGS="$LIBPATH $LIBS -lc -lnss_files -lnss_dns -lresolv" #-qstaticlink -qarch=qp

export CFLAGS="-qhot=novector -qsimd=auto -g -O3 -DUSE_MPI -DSCALAPACK $INCLUDE $DFLAGS"
export  CXXFLAGS=" -g -O3 -DUSE_MPI -DSCALAPACK $INCLUDE $DFLAGS"

./configure  --with-xerces-prefix=/nas/longleaf/home/dyost/software/qbox/xerces-c-src_2_8_0 --with-lapack=/nas/longleaf/apps/intel/17.2/intel/compilers_and_libraries_2017.2.174/linux/mkl/lib/intel64/libmkl_lapack95_lp64.a --prefix=/nas/longleaf/home/dyost/software/qball/qball_DCY/
