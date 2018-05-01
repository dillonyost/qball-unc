/* config.h.  Generated from config.h.in by configure.  */
/* config.h.in.  Generated from configure.ac by autoheader.  */

/* Define if building universal (internal helper macro) */
/* #undef AC_APPLE_UNIVERSAL_BUILD */

/* architecture */
#define ARCH "autotools"

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
/* #undef FC_DUMMY_MAIN */

/* Define if F77 and FC dummy `main' functions are identical. */
/* #undef FC_DUMMY_MAIN_EQ_F77 */

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
#define FC_FUNC(name,NAME) name ## _

/* As FC_FUNC, but for C identifiers containing underscores. */
#define FC_FUNC_(name,NAME) name ## _

/* no FFT library */
/* #undef FFT_NOLIB */

/* Define to 1 if you have the <bgpm/include/bgpm.h> header file. */
/* #undef HAVE_BGPM_INCLUDE_BGPM_H */

/* whether BlueGene/Q libraries are available */
/* #undef HAVE_BGQLIBS */

/* Defined if you have BLACS library. */
#define HAVE_BLACS 1

/* Define if you have a BLAS library. */
#define HAVE_BLAS 1

/* Define if you have a ESSL_FFT library. */
/* #undef HAVE_ESSL_FFT */

/* Define if you have a FFTW2 library. */
/* #undef HAVE_FFTW2 */

/* Define if you have a FFTW3 library. */
#define HAVE_FFTW3 1

/* Define to 1 if fseeko (and presumably ftello) exists and is declared. */
#define HAVE_FSEEKO 1

/* Define to 1 if you have the `gettimeofday' function. */
#define HAVE_GETTIMEOFDAY 1

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Defined if you have LAPACK library. */
#define HAVE_LAPACK 1

/* Define to 1 if you have the `mass' library (-lmass). */
/* #undef HAVE_LIBMASS */

/* Define if you have a MASSV library. */
/* #undef HAVE_MASSV */

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if you have the `memset' function. */
#define HAVE_MEMSET 1

/* Define to 1 if you have the `mkdir' function. */
#define HAVE_MKDIR 1

/* Define to 1 if you have the <mpi.h> header file. */
#define HAVE_MPI_H 1

/* Define to 1 if you have the `MPI_Init' function. */
#define HAVE_MPI_INIT 1

/* Define if OpenMP is enabled */
#define HAVE_OPENMP 1

/* Define to 1 if you have the `pow' function. */
#define HAVE_POW 1

/* Define if you have POSIX threads libraries and header files. */
#define HAVE_PTHREAD 1

/* Have PTHREAD_PRIO_INHERIT. */
#define HAVE_PTHREAD_PRIO_INHERIT 1

/* Defined if you have SCALAPACK library version 2 or greater. */
#define HAVE_SCALAPACK 1

/* Define to 1 if you have the <spi/include/kernel/location.h> header file. */
/* #undef HAVE_SPI_INCLUDE_KERNEL_LOCATION_H */

/* Define to 1 if you have the `sqrt' function. */
#define HAVE_SQRT 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the `strspn' function. */
#define HAVE_STRSPN 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/time.h> header file. */
#define HAVE_SYS_TIME_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the `uname' function. */
#define HAVE_UNAME 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* have the Xerces-C library */
#define HAVE_XERCES 1

/* Name of package */
#define PACKAGE "qball"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "xavier@llnl.gov"

/* Define to the full name of this package. */
#define PACKAGE_NAME "qball"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "qball alsos"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "qball"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "alsos"

/* Define to necessary symbol if this constant uses a non-standard name on
   your system. */
/* #undef PTHREAD_CREATE_JOINABLE */

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* enable MPI support */
#define USE_MPI 1

/* Version number of package */
#define VERSION "alsos"

/* the architecture is bigendian */
/* #undef WORDS_BIGENDIAN */

/* Enable large inode numbers on Mac OS X 10.5.  */
#ifndef _DARWIN_USE_64_BIT_INODE
# define _DARWIN_USE_64_BIT_INODE 1
#endif

/* Number of bits in a file offset, on hosts where this is settable. */
/* #undef _FILE_OFFSET_BITS */

/* Define to 1 to make fseeko visible on some hosts (e.g. glibc 2.2). */
/* #undef _LARGEFILE_SOURCE */

/* Define for large files, on AIX-style hosts. */
/* #undef _LARGE_FILES */

/* Define for Solaris 2.5.1 so the uint64_t typedef from <sys/synch.h>,
   <pthread.h>, or <semaphore.h> is not used. If the typedef were allowed, the
   #define below would cause a syntax error. */
/* #undef _UINT64_T */

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
/* #undef inline */
#endif

/* Define to `long int' if <sys/types.h> does not define. */
/* #undef off_t */

/* Define to `unsigned int' if <sys/types.h> does not define. */
/* #undef size_t */

/* Define to the type of an unsigned integer type of width exactly 64 bits if
   such a type exists and the standard includes do not define it. */
/* #undef uint64_t */
