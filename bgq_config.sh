#!/bin/bash

LIBHOME=/soft/libraries/alcf/current/xl
BLASDIR=$LIBHOME/BLAS/lib
LAPACKDIR=$LIBHOME/LAPACK/lib
SCALAPACK_DIR=$LIBHOME/SCALAPACK
SCALAPACKLIB=$SCALAPACK_DIR/lib/libscalapack.a
LAPACKLIB=$LAPACKDIR/liblapack.a
ESSLDIR=/soft/libraries/essl/current
HPMLIBS="-L/soft/perftools/hpctw/lib -lmpihpm_smp -L/bgsys/drivers/ppcfloor/bgpm/lib -lbgpm"
XERCESCDIR=/home/dyost/Dsource/Dlib/xerces-c-3.1.1/Dbuild/include
XERCESCLIBDIR=/home/dyost/Dsource/Dlib/xerces-c-3.1.1/Dbuild/lib
XERCESLIB=$XERCESCLIBDIR/libxerces-c.a

BGQ_SDK_PATH=/bgsys/drivers/ppcfloor
export CXX=$BGQ_SDK_PATH/comm/bin/xl/mpic++
export CC=$BGQ_SDK_PATH/comm/bin/xl/mpic++
export LIBS_BLAS="-L/soft/libraries/essl/current/5.1/lib64/ -lesslsmpbg -lesslbg"

DFLAGS=" -DPRINTALL -DUSE_ESSL -DUSE_CSTDIO_LFS -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DHPM -DUSE_XERCES -DXERCESC_3"

#INCLUDE=" -I$ESSLDIR/include -I$XERCESCDIR"

#LIBPATH=" -L$LAPACKDIR -L$BLASDIR -L$ESSLDIR/lib64 -L$XERCESCLIBDIR -L/soft/compilers/ibmcmp-aug2015/xlsmp/bg/3.1/bglib64 -L/soft/compilers/ibmcmp-aug2015/xlf/bg/14.1/bglib64"

#export LIBS="$SCALAPACKLIB $HPMLIBS -lesslsmpbg -llapack -lxlf90_r -lxlsmp -lxlfmath -lxerces-c"
#export LIBS="$SCALAPACKLIB $HPMLIBS $LAPACKLIB -lxerces-c -lxlfmath" #-lxlf90_r -lxlfmath -lxlsmp"
#export LDFLAGS="$LIBPATH $LIBS -qarch=qp -qstaticlink -lc -lnss_files -lnss_dns -lresolv"

export CFLAGS="-qhot=novector -qsimd=auto -g -O3 -DUSE_MPI -DSCALAPACK $INCLUDE $DFLAGS"
export  CXXFLAGS=" -g -O3 -qarch=qp -qsmp=omp -qtune=qp -DUSE_MPI -DSCALAPACK $INCLUDE $DFLAGS"

#./configure --with-essl-prefix=/soft/libraries/essl/current/essl/5.1/lib64/ --with-xerces-prefix=/home/dyost/Dsource/Dlib/xerces-c-3.1.1/Dbuild/ --with-lapack=/soft/libraries/alcf/current/gcc/LAPACK/lib/liblapack.a --with-scalapack=/soft/libraries/alcf/current/gcc/SCALAPACK/lib/libscalapack.a --prefix=/home/dyost/qball/qball_etrs_proj/
