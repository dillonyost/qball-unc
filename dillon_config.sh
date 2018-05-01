#!/bin/bash
# configure qball for ALCF #

LIBHOME=/nas/longleaf/apps/intel/17.2/intel/compilers_and_libraries_2017.2.174/linux/mkl/lib/intel64
MKLDIR=/nas/longleaf/apps/intel/17.2/intel/compilers_and_libraries_2017.2.174/linux/mkl/
BLASDIR=$LIBHOME
LAPACKDIR=$LIBHOME
SCALAPACK_DIR=$LIBHOME
SCALAPACKLIB=$SCALAPACK_DIR/libmkl_scalapack_ilp64.a
LAPACKLIB=$LAPACKDIR/libmkl_lapack95_ilp64.a
#ESSLDIR=/soft/libraries/essl/current
#HPMLIBS="-L/soft/perftools/hpctw/lib -lmpihpm_smp -L/bgsys/drivers/ppcfloor/bgpm/lib -lbgpm"
XERCESCDIR=/nas/longleaf/home/dyost/software/qbox/xerces-c-src_2_8_0/include
XERCESCLIBDIR=/nas/longleaf/home/dyost/software/qbox/xerces-c-src_2_8_0/lib
XERCESLIB=$XERCESCLIBDIR/libxerces-c.so

#BGQ_SDK_PATH=/bgsys/drivers/ppcfloor
export CXX=icc
export CC=mpicxx
export LIBS_BLAS="-L$LIBHOME/libmkl_blas95_ilp64.a "

DFLAGS=" -DPRINTALL -DUSE_CSTDIO_LFS -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DHPM -DUSE_XERCES -DXERCESC_2"

INCLUDE=" -I$MKLDIR/include -I$XERCESCDIR"

LIBPATH=" -L$LAPACKDIR -L$BLASDIR -L$MKLDIR/lib/intel64 -L$XERCESCLIBDIR "

export LIBS="$SCALAPACKLIB -lmkl_intel_lp64 -lmkl_lapack_lp64 -mkl_sequential -lmkl_core -lirc -lifcore -lsvml -lxerces-c  -llapack -lpthread "
#export LIBS="$SCALAPACKLIB $HPMLIBS $LAPACKLIB -lxerces-c -lxlfmath" #-lxlf90_r -lxlfmath -lxlsmp"
export LDFLAGS="$LIBPATH $LIBS -lc -lnss_files -lnss_dns -lresolv" #-qstaticlink -qarch=qp

export CFLAGS="-qhot=novector -qsimd=auto -g -O3 -DUSE_MPI -DSCALAPACK $INCLUDE $DFLAGS"
export  CXXFLAGS=" -g -O3 -DUSE_MPI -DSCALAPACK $INCLUDE $DFLAGS"

./configure  --with-xerces-prefix=/nas/longleaf/home/dyost/software/qbox/xerces-c-src_2_8_0 --with-lapack=/nas/longleaf/apps/intel/17.2/intel/compilers_and_libraries_2017.2.174/linux/mkl/lib/intel64/libimkl_lapack95_ilp64.a --with-scalapack=/nas/longleaf/apps/intel/17.2/intel/compilers_and_libraries_2017.2.174/linux/mkl/lib/intel64/libmkl_scalapack_ilp64.a --prefix=/nas/longleaf/home/dyost/qball/qball_efield/ --with-blacs=/nas/longleaf/apps/intel/17.2/intel/compilers_and_libraries_2017.2.174/linux/mkl/lib/intel64/libmkl_scalapack_ilp64.a
