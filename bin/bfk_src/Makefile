#makefile for program bfk
#  Operating system
OSTYPE=$(shell uname)
#
FC=gfortran
LINKER=gfortran
ifdef cross64
OSTYPE=MINGW
CC=x86_64-w64-mingw32-gcc -pipe $(ARCHFLAG)
CCC=x86_64-w64-mingw32-g++ -pipe $(ARCHFLAG)
FC=x86_64-w64-mingw32-gfortran -pipe $(ARCHFLAG)
LINKER=x86_64-w64-mingw32-gfortran -m64
endif

 # (linux windows)
ifneq (,$(findstring $(OSTYPE),Linux))
 LINUX = linux
else
ifneq (,$(findstring $(OSTYPE),Darwin))
 DARWIN = darwin
else
 WINDOWS = windows
ifneq (,$(findstring CYGWIN,$(OSTYPE)))
 CYGWIN=1
endif
endif
endif

ifdef WINDOWS
EXE=.exe
STATIC=-static
endif

all :  clean bfk.o 
#../chibfk$(EXE) 

bfk.o: bfk.f95
	$(FC) $(STATIC) -lgfortran -std=legacy -o ../bfk $^ 

#../chibfk$(EXE): chibfk.f95
#	gfortran -static -lgfortran -o ../chibfk $^ 

clean:
	rm -vf ./*.o commondata.mod matrixelements.mod
