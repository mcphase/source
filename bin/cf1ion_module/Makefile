#makefile for program so1ion

vpath %.h  ../include
vpath %.hpp ../include
vpath %.o  ../src/vector

###########################################################################
# gcc version egcs-2.95.2 19991024 (release) on a Linux system
# Intel Pentium PC, use -m486 for 486 CPU or -mpentium for Pentium Systems
###########################################################################

include ../src/Makefile.common

###########################################################################
### include directory: don't change !					  #
###########################################################################

INCLUDE= ../../include

###########################################################################
### general rules							  #
###########################################################################

.SUFFIXES: .cc .o

#	${CC} -c $(patsubst -Wall,,${COPTS}) -I ${INCLUDE} $<
.c.o: 
	${CC} -c ${COPTS} -I ${INCLUDE} $<

.cc.o: 
	${CCC} -c ${CCOPTS}   -I ${INCLUDE} $<

###########################################################################
###########################################################################



INCLDIRS =      -I. \
		-I../include \
		-I/usr/include/readline \
		-I/local/lib  \
		-I/usr/X11R6/include

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

all : cfieldrout.o cfield.so ../cfield$(EXE) ../so1ion$(EXE) module.o

VECOBJS= matpack.o sortdbl.o sort2dbl.o sort3dbl.o ivector.o \
		imatrix.o dvector.o dmatrix.o dmatmul.o cvector.o \
		cmatrix.o cmatmul.o dludcomp.o dsvdcomp.o dortho.o \
		cludcomp.o tred.o imtql.o imtql2.o rseigen.o cheigen.o \
		chgeigen.o matlapl.o moment.o matexph.o matsqth.o magicsquare.o

CFOBJS=  diahermx.o eingabe.o intensit.o komplex.o matrix.o mb11a.o minima.o \
	       ortho.o spline.o stevens.o theta.o va05a.o


${CFOBJS}: types.c
cfield.o: types.c main.c
so1ion.o: types.c main.c

# used for internal module in mcphas and mcdisp
cfieldrout.o: cfieldrout.c $(CFOBJS) 

cfield.so: module.o cfield.o  $(VECOBJS) $(CFOBJS) 
	$(CCC) $(CCOPTS) -Wall -o cfield.so $^ -shared $(LFLAGS)

../cfield$(EXE): cfield.o $(CFOBJS)
	$(CCC) $(CCOPTS) -Wall -o ../cfield $^ $(LFLAGS)

../so1ion$(EXE): so1ion.o $(CFOBJS)
	$(CCC) $(CCOPTS) -Wall -o ../so1ion $^ $(LFLAGS)

module.o: module.c
	$(CCC) $(CCOPTS) -Wall -fPIC -c  -o module.o module.c  $(INCLDIRS) $(LFLAGS)

clean:
	rm -vf ./*.o ./*.info ./*.parameter ./*.output

