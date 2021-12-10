#*********************************************************************
#
# File:          Makefile
#
# Project:       McPhase - CF Phonon Module bcfph
#
# Description:   Makefile for bcfph  - compiles and links all c++ code
#
# Author(s):     M. Rotter
#
# Original:	 February 2021
#
#**********************************************************************

# make looks in the following directories for files
vpath %.c  ./
vpath %.cpp ./
vpath %.h  ../../bin/include ../../bin/ic1ion_module
vpath %.hpp ../../bin/include ../../bin/ic1ion_module
vpath %.o  ./ ../../bin/src ../../bin/src/vector ../../bin/cf1ion_module ../../bin/src/functions

#---------------------------------------------------------------------
#  Operating system
OSTYPE=$(shell uname)
#  (linux windows)
ifneq (,$(findstring $(OSTYPE),Linux))
 LINUX = linux
else
ifneq (,$(findstring $(OSTYPE),Darwin))
 LINUX = linux
else
 WINDOWS = windows
ifneq (,$(findstring CYGWIN,$(shell uname)))
 CYGWIN=1
endif
endif
endif
#
 STATIC = static
#
#---------------------------------------------------------------------
#               ... and his Flags
#
#                Includes
#

INCLDIRS =      -I. \
		-I../../bin/include \
		-I../../bin/ic1ion_module \
#		-I/usr/include/readline \
		-I/local/lib  \
		-I/usr/X11R6/include

ifdef CYGWIN
CYGFLAGS = -DCYGWIN
endif 

include ../../bin/src/Makefile.common

CFLAGS  = -c $(CCOPTS) $(INCLDIRS) $(CYGFLAGS) -U_FORTIFY_SOURCE

# PGPLLIBS =       -L$(PGPLOT_DIR)  -lcpgplot -lpgplot
ifdef LINUX
XLIBS= -ldl
#XLIBS=           -L/usr/X11/lib -lX11 -lXext -lc -lXpm -ldl
#EZLIB=           -ljpeg -ldl
else
EXE=.exe
endif
ifdef OMP
LFLAGS +=          $(XLIBS) $(EZLIB) -lstdc++ -lm -lgomp -lpthread
else
#LFLAGS =          $(XLIBS) $(EZLIB) -lstdc++ -lm
LFLAGS +=          $(XLIBS) $(EZLIB) -lstdc++ -lm ${LBLAS}
endif

# uncomment next line and set NUM_THREADS in mcphas_htcalc.c and mcdisp.c
# according to your machine if you want to use multiple processing to speed up
# calculations
threads=1
# 
#ifdef threads
THREADFLAG=-D_THREADS
ifdef LINUX
THREADLINK=-lpthread
endif
#endif
#---------------------------------------------------------------------
# RCS options
#
#               RCS lock options
RCSLOCK  =      co -l -r$(VERSION)
#               RCS check out options
RCSCO    =      co -r$(VERSION)
#               RCS check in options
RCSCI    =      ci -u -f -s"Rel" -r$(VERSION) -m"$(LOCKMSG)"


#---------------------------------------------------------------------
# All files needed for the library
#
#	File list for RCS control
#               all include files

#               source files


#--------------------------------------------------------------------
# What has to be made
#
#               Names of the include files in the 
#		installation directories
#

INCL = mcphas.h  testspin.hpp qvectors.hpp spincf.hpp par.hpp  jjjpar.hpp inipar.hpp \
       mfcf.hpp vector.h martin.h graphic_parameters.hpp
OBJS = jjjpar.o par.o physprop.o qvectors.o testspin.o spincf.o inipar.o \
        mfcf.o martin.o ionpars.o graphic_parameters.o perlparse.o

VOBS = ivector.o imatrix.o dvector.o dmatrix.o dmatmul.o cvector.o cmatrix.o cmatmul.o \
	dludcomp.o dsvdcomp.o dortho.o cludcomp.o tred.o imtql.o imtql2.o rseigen.o cheigen.o \
	chgeigen.o matlapl.o moment.o matexph.o matsqth.o magicsquare.o \
	matpack.o sortdbl.o sort2dbl.o sort3dbl.o


#---------------------------------------------------------------------
# Build bcfph

all :   ../bcfph$(EXE)

# ****************************************************
# compile sources 


bcfph.o: bcfph.c martin.h myev.h
	$(CCC) $(CFLAGS) $(THREADFLAG) $< 


# ************************************************
# make  binaries (link .o files)


../bcfph$(EXE):  bcfph.o martin.o  myev.o  $(VOBS)
	@echo Linking ../../bin/bcfph
	$(LINKER)  -Wall  -o ../../bin/bcfph $^  $(LFLAGS) $(THREADLINK)	

clean: 
	rm -vf ./*.o

