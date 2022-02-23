#makefile for program bfk


#  Operating system
OSTYPE=$(shell uname)
#  (linux windows)
ifneq (,$(findstring $(OSTYPE),Linux))
 LINUX = linux
else
ifneq (,$(findstring $(OSTYPE),Darwin))
 DARWIN = darwin
else
 WINDOWS = windows
ifneq (,$(findstring CYGWIN,$(shell uname)))
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
	gfortran $(STATIC) -lgfortran -o ../bfk $^ 

#../chibfk$(EXE): chibfk.f95
#	gfortran -static -lgfortran -o ../chibfk $^ 

clean:
	rm -vf ./*.o commondata.mod matrixelements.mod
