#makefile for program cowan
FC=gfortran
LINKER=gfortran
ifdef cross64
OSTYPE=MINGW
CC=x86_64-w64-mingw32-gcc -pipe $(ARCHFLAG)
CCC=x86_64-w64-mingw32-g++ -pipe $(ARCHFLAG)
FC=x86_64-w64-mingw32-gfortran -pipe $(ARCHFLAG)
LINKER=x86_64-w64-mingw32-gfortran -m64
endif


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
ifneq (,$(findstring CYGWIN,$(OSTYPE)))
 CYGWIN=1
endif
endif
endif

ifdef WINDOWS
EXE=.exe
STATIC=-static
endif

all :  clean  RCN36K.o RCN2K.o
# ../RCG11K$(EXE)   not needed
#../chibfk$(EXE) 

RCN2K.o: RCN2K.F
	$(FC) $(STATIC) -lgfortran -o ../RCN2K $^ 

RCN36K.o: RCN36K.F
	$(FC) $(STATIC) -lgfortran -o ../RCN36K $^ 

#../RCG11K$(EXE): RCG11K.F
#	$(LINKER) $(STATIC) -lgfortran -o ../RCG11K $^ 

#../chibfk$(EXE): chibfk.f95
#	$(LINKER) -static -lgfortran -o ../chibfk $^ 

clean:
	rm -vf ./*.o
