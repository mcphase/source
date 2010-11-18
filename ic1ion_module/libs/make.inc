####################################################################
#  LAPACK make include file.                                       #
#  LAPACK, Version 3.1.1                                           #
#  February 2007                                                   #
####################################################################
#
# See the INSTALL/ directory for more examples.
#
SHELL = /bin/sh
#
#  The machine (platform) identifier to append to the library names
#
PLAT = _LINUX
#  
#  Modify the FORTRAN and OPTS definitions to refer to the
#  compiler and desired compiler options for your machine.  NOOPT
#  refers to the compiler options desired when NO OPTIMIZATION is
#  selected.  Define LOADER and LOADOPTS to refer to the loader and 
#  desired load options for your machine.
#
FORTRAN  = gfortran
ifeq ($(shell uname),Linux)
 LINUX = LINUX
 PLAT = _LINUX
 FORTRAN=gfortran --std=legacy
 LIBFORTRAN=-lgfortran
else
 MSYS=1
 PLAT = _WINDOWS
 FORTRAN= g77
 LIBFORTRAN=-lg2c
endif

ifneq ($(shell gcc --version|grep " 3."),)
 FORTRAN=g77
 LIBFORTRAN=-lg2c
else
ifneq ($(shell gcc --version|grep " 4."),)
 FORTRAN=gfortran --std=legacy
 LIBFORTRAN=-lgfortran
endif
endif

ifneq (,$(findstring CYGWIN,$(shell uname)))
ifneq ($(shell gcc --version|grep " 3."),)
FLAGS+=-mno-cygwin
endif
ifdef cross64
FORTRAN  = x86_64-w64-mingw32-gfortran --std=legacy
ARCH     = x86_64-w64-mingw32-ar
ARCHFLAGS= cr
RANLIB   = x86_64-w64-mingw32-ranlib
endif
ifdef cross32
FORTRAN  = i686-w64-mingw32-gfortran --std=legacy
ARCH     = i686-w64-mingw32-ar
ARCHFLAGS= cr
RANLIB   = i686-w64-mingw32-ranlib
endif
endif

ifdef MSYS
FLAGS   += -pipe -g
else
FLAGS   += -pipe -g -fPIC
endif
ifdef fast
ifneq ($(shell gcc --version|grep " 3."),)
FLAGS   += -O2
else
FLAGS   += -O2 -march=native
endif
ifdef profile
FLAGS   += -pg
else
FLAGS   += -fomit-frame-pointer
endif
endif
ifdef distrib
FLAGS   += -mtune=pentium3 -march=i686
endif
ifdef b32
FLAGS   += -m32
NOOPT   += -m32
endif
OPTS     = $(FLAGS)
DRVOPTS  = $(OPTS)
NOOPT    = -fPIC -O0
LOADER   = $(FORTRAN)
LOADOPTS =
CC       = gcc
CFLAGS   = $(FLAGS)
#
# Timer for the SECOND and DSECND routines
#
# Default : SECOND and DSECND will use a call to the EXTERNAL FUNCTION ETIME
# TIMER    = EXT_ETIME
# For RS6K : SECOND and DSECND will use a call to the EXTERNAL FUNCTION ETIME_
# TIMER    = EXT_ETIME_
# For gfortran compiler: SECOND and DSECND will use a call to the INTERNAL FUNCTION ETIME
TIMER    = INT_ETIME 
# If your Fortran compiler does not provide etime (like Nag Fortran Compiler, etc...)
# SECOND and DSECND will use a call to the INTERNAL FUNCTION CPU_TIME 
# TIMER    = INT_CPU_TIME
# If neither of this works...you can use the NONE value... In that case, SECOND and DSECND will always return 0
# TIMER     = NONE
#
#  The archiver and the flag(s) to use when building archive (library)
#  If you system has no ranlib, set RANLIB = echo.
#
ifndef ARCH
ARCH     = ar
ARCHFLAGS= cr
RANLIB   = ranlib
endif
#
#  The location of the libraries to which you will link.  (The 
#  machine-specific, optimized BLAS library should be used whenever
#  possible.)
#
BLASLIB      = libblas.a
LAPACKLIB    = liblapack.a
ARPACKLIB    = libarpack.a
