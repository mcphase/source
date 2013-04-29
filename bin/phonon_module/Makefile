#makefile for dynamically loaded module

vpath %.h  ../include
vpath %.hpp ../include
vpath %.o  ../src/vector

INCLDIRS =      -I. \
		-I../include \
		-I/usr/include/readline \
		-I/local/lib  \
		-I/usr/X11R6/include


all : phonon.so

VECOBJS= matpack.o sortdbl.o sort2dbl.o sort3dbl.o ivector.o \
		imatrix.o dvector.o dmatrix.o dmatmul.o cvector.o \
		cmatrix.o cmatmul.o dludcomp.o dsvdcomp.o dortho.o \
		cludcomp.o tred.o imtql.o imtql2.o rseigen.o cheigen.o \
		chgeigen.o matlapl.o moment.o matexph.o matsqth.o magicsquare.o

phonon.so: phonon.o  $(VECOBJS)
	g++  -Wall -o phonon.so $^ -shared -static-libgcc -static-libgfortran

phonon.o: phonon.c
	g++  -Wall -fPIC -c  -o phonon.o phonon.c  $(INCLDIRS)

clean:
	rm -vf phonon.o

