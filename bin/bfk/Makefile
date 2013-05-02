#makefile for program bfk


#  Operating system
OSTYPE=$(shell uname)
#  (linux windows)
ifneq (,$(findstring $(OSTYPE),Linux))
 LINUX = linux
else
 WINDOWS = windows
ifneq (,$(findstring CYGWIN,$(shell uname)))
 CYGWIN=1
endif
endif

ifdef WINDOWS
EXE=.exe
endif

all :  clean ../bfk$(EXE) ../chibfk$(EXE) 

../bfk$(EXE): bfk.f95
	gfortran -static -lgfortran -o ../bfk $^ 

../chibfk$(EXE): chibfk.f95
	gfortran -static -lgfortran -o ../chibfk $^ 

clean:
	rm -vf ./*.o commondata.mod matrixelements.mod

        