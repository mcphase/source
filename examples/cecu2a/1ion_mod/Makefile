#makefile for dynamically loaded module
#for single ion properties of Ce in CeCu2
vpath %.h  ../../../bin/include
vpath %.hpp ../../../bin/include
vpath %.o  ../../../bin/src/vector

include ../../../bin/src/Makefile.common

INCLDIRS =      -I. \
		-I../../../bin/include \
		-I../../../bin/ic1ion_module \
		-I/usr/include/readline \
                -I/Library/Developer/CommandLineTools/SDKs/MacOSX15.2.sdk/System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/Headers/Sparse/ \
		-I/local/lib  \
		-I/usr/X11R6/include


all : kramer.so fit

VECOBJS= matpack.o sortdbl.o sort2dbl.o sort3dbl.o ivector.o \
		imatrix.o dvector.o dmatrix.o dmatmul.o cvector.o \
		cmatrix.o cmatmul.o dludcomp.o dsvdcomp.o dortho.o \
		cludcomp.o tred.o imtql.o imtql2.o rseigen.o cheigen.o \
		chgeigen.o matlapl.o moment.o matexph.o matsqth.o magicsquare.o

kramer.so: kramer.o  $(VECOBJS)
	$(CCC)  -Wall -o kramer.so $^ -shared -L../../../bin/ic1ion_module/libs $(LFLAGS)
#-larpack -llapack -lblas -L/usr/local/gfortran/lib -lgfortran
# -static-libgcc -static-libgfortran -L../../../bin/ic1ion_module/libs -llapack -lblas -lgfortran
#--out-implib 

kramer.o: kramer.c 
	$(CCC) -Wall -fPIC -c  -o kramer.o kramer.c  $(INCLDIRS)

fit:
	cp -f ../fit/mcphas.jjj.sav ../fit/mcphas.jjj

clean:
	rm -vf kramer.o 
	rm -vf ../results/*
	rm -vf ../fit/mcphas.jjj.fit

