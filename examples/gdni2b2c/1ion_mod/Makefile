#makefile for dynamically loaded module
#for single ion properties of Gd in GdNi2B2C
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

VECOBJS= matpack.o sortdbl.o sort2dbl.o sort3dbl.o ivector.o \
		imatrix.o dvector.o dmatrix.o dmatmul.o cvector.o \
		cmatrix.o cmatmul.o dludcomp.o dsvdcomp.o dortho.o \
		cludcomp.o tred.o imtql.o imtql2.o rseigen.o cheigen.o \
		chgeigen.o matlapl.o moment.o matexph.o matsqth.o magicsquare.o

all : brillouin.so bosefermi.so

brillouin.so: brillouin.o  $(VECOBJS)
	$(CCC)  -Wall -o brillouin.so $^  -shared -L../../../bin/ic1ion_module/libs $(LFLAGS)
# -static-libgcc -static-libgfortran -L../../../bin/ic1ion_module/libs -llapack -lblas -lgfortran


bosefermi.so: bosefermi.o  $(VECOBJS)
	$(CCC)  -Wall -o bosefermi.so $^  -shared -L../../../bin/ic1ion_module/libs $(LFLAGS)
# -larpack -llapack -lblas -L/usr/local/gfortran/lib -lgfortran
# -static-libgcc -static-libgfortran  -L../../../bin/ic1ion_module/libs -llapack -lblas -lgfortran


brillouin.o: brillouin.c 
	$(CCC)  -Wall -c  -fPIC -std=c++0x -o brillouin.o brillouin.c  $(INCLDIRS)

bosefermi.o: bosefermi.c 
	$(CCC)  -Wall -c  -fPIC -std=c++0x -o bosefermi.o bosefermi.c  $(INCLDIRS)
        

clean:
	rm -vf brillouin.o bosefermi.o

