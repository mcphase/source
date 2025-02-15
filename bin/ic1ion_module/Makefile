OSTYPE=$(shell uname)

ifdef icc
CC = icpc
LINKER = icpc
else
CC = g++
LINKER = g++
PEDANT = --pedantic
endif


# Determines if we have a gcc4.x.x compiler
ifeq ($(findstring false,$(shell ./checkgcc.sh -4 $(CC))),false)
GCC4=1
endif
ifeq ($(findstring false,$(shell ./checkgcc.sh -4.2 $(CC))),false)
GCC42=1
endif
ifeq ($(findstring false,$(shell ./checkgcc.sh -5 $(CC))),false)
GCC4=1
GCC42=1
GCC5=1
endif

ifdef cross64
OSTYPE=MINGW
CC=x86_64-w64-mingw32-g++ -m64 $(ARCHFLAG)
LINKER=x86_64-w64-mingw32-g++ -m64
endif


# Figures operating system
ifneq (,$(findstring MINGW,$(OSTYPE)))
MSYS=1
WINDOWS=1
EXE=.exe
else
ifeq ($(shell uname),Linux)
LINUX = LINUX
LIBFORTRAN=-lgfortran
else
ifneq (,$(findstring CYGWIN,$(OSTYPE)))
CYGWIN=1
WINDOWS=1
EXE=.exe
endif                # if CYGWIN
endif                # if LINUX
endif                # if MSYS

ifdef CYGWIN
# To fix an error in spincf where min and max is used as variable names
#   and conflicts with definitions in <w32api/windef.h>
CFLAGS+=-DNOMINMAX
ifndef GCC4
CFLAGS+=-mno-cygwin
LINKER+=-mno-cygwin
endif
# Checks to see if we want to use the 64bit crosscompiler in Cygwin
ifdef cross64
CC=x86_64-w64-mingw32-g++ -m64 $(ARCHFLAG)
LINKER=x86_64-w64-mingw32-g++ -m64
endif
ifdef cross32
CC=i686-w64-mingw32-g++ -m32 $(ARCHFLAG)
LINKER=i686-w64-mingw32-g++ -m32
endif
endif

ifdef GCC4
 LIBFORTRAN=-lgfortran
 LDFLAGS=-L/usr/local/gfortran/lib/ 
# LDFLAGS=-L/opt/homebrew/Cellar/gcc/11.3.0_1/lib/gcc/11/
else
 LIBFORTRAN=-lg2c
 CFLAGS+=-D_G77
endif

# Determines if we should use ARPACK or not
ifdef noarpack
CFLAGS+=-DNO_ARPACK
else
LIBARPACK=-larpack
endif

# Determines if we have GOTOBLAS
ifneq ($(strip $(wildcard libs/GotoBLAS/libgoto.a)),)
ifndef cross64
gotoblas=1
endif
endif

CFLAGS+= -fpermissive
##################################################################################
#                              Declarations                                      #
##################################################################################
#------------------------------    MSYS    --------------------------------------#
ifdef WINDOWS
CFLAGS+= -g -Wall $(PEDANT) -D_WINDOWS
else
CFLAGS+= -g -ansi -Wall $(PEDANT) -fPIC
endif
#------------------------------    FAST    --------------------------------------#
ifdef fast
ifdef icc
CFLAGS+= -xHOST -O2 -IPF_fma -long_double -IPF_fltacc#-ipo
else
ifdef GCC42
CFLAGS+= -pipe -O3 -march=native
#CFLAGS+= -O2 -pipe -fomit-frame-pointer -march=core2 -msse4 -msse4.1 -msse4.2 -mcx16 -msahf
else
CFLAGS+= -pipe -O2
endif
endif
ifndef profile
CFLAGS+= -fomit-frame-pointer
endif
endif
#------------------------------   LINKER   --------------------------------------#
ifdef matlab
ifdef MSYS
LDFLAGS = -L. -lmwrefblas -lmwlapack -lmwarpack
else
LDFLAGS = -L/usr/local/MATLAB/bin/glnxa64 -Wl,-rpath,/usr/local/MATLAB/bin/glnxa64 -lut -lmwbinder -lmwblas -lmwlapack -Llibs $(LIBARPACK) $(LIBFORTRAN)
endif
else    # if matlab
ifdef gotoblas
LDFLAGS = -Llibs -Llibs/GotoBLAS $(LIBARPACK) -llapack -lgoto $(LIBFORTRAN)
ifeq (,$(findstring CYGWIN,$(OSTYPE)))
LDFLAGS += -lpthread
endif
else    # if gotoblas
ifdef MKLROOT
LDFLAGS+= -L../ic1ion_module/libs -larpack $(LIBFORTRAN) -fopenmp -m64 -I$(MKLROOT)/include -L$(MKLROOT)/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm
CFLAGS+= -D_G77
else    # if mkl
#LDFLAGS = $(LIBARPACK) -llapack -lf77blas -latlas
ifdef cross64
LDFLAGS =  -Llibs $(LIBARPACK) -llapack -lblas  $(LIBFORTRAN) -lquadmath
else
LDFLAGS =  -Llibs $(LIBARPACK) -llapack -lblas -L/usr/local/gfortran/lib $(GFORTRANLIB) $(LIBFORTRAN)
à LDFLAGS= -Llibs $(LIBARPACK) -llapack -lblas -L/opt/homebrew/Cellar/gcc/11.3.0_1/lib/gcc/11/  $(LIBFORTRAN)
endif
#LDFLAGS = -Llibs $(LIBARPACK) -llapack -L/home/eqr/src/ATLAS-3.8.3/cygwin32_penryn/lib -lf77blas -latlas -lgfortran
#LDFLAGS = -Llibs $(LIBARPACK) -llapack -L/home/eqr/src/ATLAS-3.8.3/mingw64_penryn/lib -lf77blas -latlas -lgfortran
ifdef MSYS
ifdef GCC5
LDFLAGS+= -lquadmath
endif   
ifdef GCC4
LDFLAGS+= -static
endif   
endif   # ifdef MSYS
endif   # else ~mkl
endif   # else ~gotoblas
endif   # else ~matlab
ifdef icc
LDFLAGS += -lifcore
endif
#------------------------------  32/64bit  --------------------------------------#
ifdef b32
LINKER += -m32
CFLAGS += -m32
LDFLAGS += -m32
endif
#----------------------------  GNU Profiler  ------------------------------------#
ifdef profile
CFLAGS+= -pg
LDFLAGS+= -pg
endif
#------------------------------   DISTRO   --------------------------------------#
ifdef distrib
#CFLAGS += -mtune=pentium3 -march=i686
STATIC = -static
endif
##################################################################################

# List of object files to be created and their headers
OBJECT = maths               icpars icmf par_io ic_hmltn ic_diag mmio njsyms cfp coulomb so_cf lovesey spectre
HEADER = maths states ic1ion icpars icmf

.PHONY: default all clean

default: all
all: ../ic1ionit$(EXE) ic1ion.so ../icf1ionit$(EXE) icf1ion.so 

# Defines some variables
OBJ := $(foreach objs,$(OBJECT),$(objs).o)
HPP := $(foreach head,$(HEADER),$(head).hpp)

VECOBJ = ../src/vector/cheigen.o ../src/vector/chgeigen.o ../src/vector/cludcomp.o \
         ../src/vector/cmatmul.o ../src/vector/cmatrix.o ../src/vector/cvector.o \
         ../src/vector/dludcomp.o ../src/vector/dmatmul.o ../src/vector/dmatrix.o \
         ../src/vector/dortho.o ../src/vector/dsvdcomp.o ../src/vector/dvector.o \
         ../src/vector/imatrix.o ../src/vector/imtql.o ../src/vector/imtql2.o \
         ../src/vector/ivector.o ../src/vector/magicsquare.o ../src/vector/matexph.o \
         ../src/vector/matlapl.o ../src/vector/matpack.o ../src/vector/matsqth.o \
         ../src/vector/moment.o ../src/vector/rseigen.o ../src/vector/sort2dbl.o \
         ../src/vector/sort3dbl.o ../src/vector/sortdbl.o ../src/vector/tred.o

#ifeq ($(findstring $(shell ./checkgcc.sh -4 $(CC)),"false"),)
FORCE:

checkgcc: FORCE
	./checkgcc.sh 3.4.5 $(CC) -warn

libs/liblapack.a:
	cd libs && make

ifdef integral
../ic1ionit$(EXE): checkgcc libs/liblapack.a $(OBJ) states.o ic1ion.o truncate.o ic1ion_module.o  ../src/martin.o
	@echo Linking $@
	$(LINKER) -o ../ic1ion ic1ion.o truncate.o states.o ic1ion_module.o  ../src/martin.o $(OBJ) $(VECOBJ) $(LDFLAGS) $(STATIC)
else
../ic1ionit$(EXE): checkgcc libs/liblapack.a $(OBJ) states.o ic1ion.o  ../src/martin.o
	@echo Linking $@
	$(LINKER) -o ../ic1ion ic1ion.o states.o  ../src/martin.o $(OBJ) $(VECOBJ) $(LDFLAGS) $(STATIC)
endif

ic1ion.o: ic1ion.cpp maths.hpp states.hpp ic1ion.hpp
ifdef integral
	$(CC) $(CFLAGS) -I../include -I. -D_INTEGRAL -c $<
else
	$(CC) $(CFLAGS) -I../include -I. -c $<
endif

states.o: states.cpp states.hpp ../include/constants.h
	$(CC) $(patsubst -O3,   ,${CFLAGS}) -I../include -c $<

truncate.o: truncate.cpp
	$(CC) $(CFLAGS) -I../include -I. -c $<

ic1ion_module.o: ic1ion_module.cpp ic1ion_module.hpp ../include/singleion_module.hpp
ifdef integral
	$(CC) $(CFLAGS) -D_INTEGRAL -I../include -I. -c $<
else
ifdef WINDOWS
	$(CC) $(CFLAGS) -I../include -I. -DBUILD_DLL -c $<
else
	$(CC) $(CFLAGS) -I../include -I. -c $<
endif
endif


icf1ion_module.o: icf1ion_module.cpp icf1ion_module.hpp ../include/sparsecomplex.hpp
ifdef integral
	$(CC) $(CFLAGS) -D_INTEGRAL -I../include -I. -c $<
else
	$(CC) $(CFLAGS) -I../include -I. -c $<
endif

ic1ion.so: checkgcc libs/liblapack.a $(OBJ) states.o ic1ion.o truncate.o ic1ion_module.o  ../src/martin.o
	@echo Linking $@
ifdef WINDOWS
	@$(LINKER) -o ic1ion.so ic1ion_module.o ic1ion.o truncate.o states.o  ../src/martin.o $(OBJ) $(VECOBJ) -shared $(LDFLAGS) -Wl,--out-implib,libic1ion_module.a
else
	@$(LINKER) -o ic1ion.so ic1ion_module.o ic1ion.o truncate.o states.o  ../src/martin.o $(OBJ) $(VECOBJ) -shared $(LDFLAGS)
endif



ICFOBJ := par_io.o icpars.o cfp.o ic_diag.o maths.o njsyms.o states.o mmio.o ../src/martin.o # lovesey.o icmf.o ic_hmltn.o coulomb.o so_cf.o
../icf1ionit$(EXE): checkgcc libs/liblapack.a $(ICFOBJ) icf1ion_module.o
	@echo Linking $@
	@$(LINKER) -o ../icf1ion icf1ion_module.o $(ICFOBJ) $(VECOBJ) $(LDFLAGS) $(STATIC)

icf1ion.so: checkgcc libs/liblapack.a $(ICOBJ) icf1ion_module.o
	@echo Linking $@
ifdef WINDOWS
	@$(LINKER) -o icf1ion.so icf1ion_module.o $(ICFOBJ) $(VECOBJ) -shared $(LDFLAGS) -Wl,--out-implib,libicf1ion_module.a
else
	@$(LINKER) -o icf1ion.so icf1ion_module.o $(ICFOBJ) $(VECOBJ) -shared $(LDFLAGS)
endif

icmf.o : icmf.cpp ic1ion.hpp ../include/constants.h 
	$(CC) $(CFLAGS) -I../include -I. -c $<

# Pattern rules
# Syntax is
# targets ...: target-pattern: dep-patterns ...
#	commands
#	...

$(OBJ): %.o: %.cpp $(HPP) 
	$(CC) $(CFLAGS) -I../include -I. -c $<

# Implicit rules (suffix rules)
# NB. $@ = current target, $< = implied source, $* = filename without suffix.
# Syntax is
# .insuffix.outsuffix:
#	command
# For the implicit rules to work we must declare the suffixes used here
#.SUFFIXES: .tex .aux .toc .bbl .blg .log .dvi .ps .pdf .myt
#.tex.pdf:
#ifneq ($(shell if [ $*.aux -nt $*.bbl ]; then echo 1; fi),1)
#	$(LaTeX) $< $(endnot)
#endif
#	$(LaTeX) $< $(endnot)
#
#.aux.bbl:
#	$(BibTeX)  $* $(endnot)
#
#.tex.aux:
#	$(LaTeX) $< $(endnot)

clean:
	rm -vf *.o
	rm -vf *~
	rm -vf checkgcc checkgcc.exe

cleanall: clean
	cd libs && make clean
