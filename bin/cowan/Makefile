#makefile for program cowan


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

all :  clean  ../RCN36K$(EXE) ../RCN2K$(EXE)
# ../RCG11K$(EXE)   not needed
#../chibfk$(EXE) 

 ../RCN2K$(EXE): RCN2K.F
	gfortran $(STATIC) -lgfortran -o ../RCN2K $^ 

../RCN36K$(EXE): RCN36K.F
	gfortran $(STATIC) -lgfortran -o ../RCN36K $^ 

#../RCG11K$(EXE): RCG11K.F
#	gfortran $(STATIC) -lgfortran -o ../RCG11K $^ 

#../chibfk$(EXE): chibfk.f95
#	gfortran -static -lgfortran -o ../chibfk $^ 

clean:
	rm -vf ./*.o
