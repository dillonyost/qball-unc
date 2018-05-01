#-------------------------------------------------------------------------------
#
#  BlueGene/Q (IBM Q32)
#
#-------------------------------------------------------------------------------
# $Id: linux_x86_64_intel.mk,v 1.1 2010/01/22 00:35:01 draeger1 Exp $
#
 PLT=BGQ
#-------------------------------------------------------------------------------

# do this to force manual compilation of key threaded objects
 OMPHACK = 1

 LIBHOME =/soft/libraries/alcf/current/xl
 BLASDIR=$(LIBHOME)/BLAS/lib
 LAPACKDIR=$(LIBHOME)/LAPACK/lib
 SCALAPACK_DIR = $(LIBHOME)/SCALAPACK
 SCALAPACKLIB  = $(SCALAPACK_DIR)/lib/libscalapack.a
 # DCY ____ ESSLDIR =/soft/libraries/essl/current/essl/5.1
 ESSLDIR = /soft/libraries/essl/current/
 HPMLIBS = -L/soft/perftools/hpctw/lib -lmpihpm_smp -L/bgsys/drivers/ppcfloor/bgpm/lib -lbgpm
 JAGGEMMLIB =# $(LIBHOME)/jaggemm_opt/libjaggemm.a
 #CTFDIR = $(LIBHOME)/ctf-latest/cyclopstf
 #CTFLIB = -L$(LIBHOME)/lib -lcyclopstf.jag
 XERCESCDIR=/home/dyost/Dsource/Dlib/xerces-c-3.1.1/Dbuild/include
 XERCESCLIBDIR=/home/dyost/Dsource/Dlib/xerces-c-3.1.1/Dbuild/lib
 XERCESLIB=$(XERCESCLIBDIR)/libxerces-c.a

 BGQ_SDK_PATH =/bgsys/drivers/ppcfloor
 CXX=$(BGQ_SDK_PATH)/comm/bin/xl/mpic++
 CC=$(BGQ_SDK_PATH)/comm/bin/xl/mpic++
 #FORT=$(BGQ_SDK_PATH)/comm/gcc/bin/mpixlf90_r
 AR=#$(BGQ_SDK_PATH)/gnu-linux/lib64/bin/ar
 RANLIB=#$(BGQ_SDK_PATH)/gnu-linux/lib64/

 LD=$(CXX)

 DFLAGS += -DPRINTALL -DUSE_ESSL -DUSE_CSTDIO_LFS \
	-D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DHPM -DUSE_XERCES -DXERCESC_3
 
 INCLUDE = -I$(ESSLDIR)/include -I$(XERCESCDIR)
 
 CXXFLAGS= -g -O3 -qarch=qp -qsmp=omp -qtune=qp -DUSE_MPI -DSCALAPACK -D$(PLT) $(INCLUDE) $(DFLAGS)
 CFLAGS= -qhot=novector -qsimd=auto -g -O3 -DUSE_MPI -DSCALAPACK -D$(PLT) $(INCLUDE) $(DFLAGS)

 LIBPATH = -L$(LAPACKDIR) -L$(BLASDIR) -L$(ESSLDIR)/lib64 -L/soft/compilers/ibmcmp-nov2014/xlsmp/bg/3.1/bglib64 \
	-L/soft/compilers/ibmcmp-nov2014/xlf/bg/14.1/bglib64 -L$(XERCESCLIBDIR)
 LIBS =  $(SCALAPACKLIB) $(HPMLIBS) -lesslsmpbg -llapack -lxlf90_r -lxlsmp -lxlfmath -lxerces-c
 LDFLAGS = $(LIBPATH) $(LIBS) -qarch=qp -qstaticlink -lc -lnss_files -lnss_dns -lresolv

#TAUROOTDIR = $(LIBHOME)/tau/tau-2.21.2
#ifneq (,$(findstring DTAU,$(DFLAGS)))
#        include  $(TAUROOTDIR)/include/Makefile
#        CXXFLAGS+=$(TAU_INCLUDE) $(TAU_DEFS)
#       LIBS+=$(TAU_MPI_LIBS) $(TAU_LIBS)                                                                      #            
#        LDFLAGS+= $(TAU_LIBS)
#endif


#-------------------------------------------------------------------------------
