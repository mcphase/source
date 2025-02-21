#makefile for dynamically loaded module
vpath %.h  $(MCPHASE_DIR)/bin/include
vpath %.hpp $(MCPHASE_DIR)/bin/include
vpath %.o  $(MCPHASE_DIR)/bin/src/vector
vpath %.o  $(MCPHASE_DIR)/bin/src

include $(MCPHASE_DIR)/bin/src/Makefile.common

INCLDIRS =      -I. \
		-I$(MCPHASE_DIR)/bin/include \
		-I$(MCPHASE_DIR)/bin/ic1ion_module \
		-I/usr/include/readline \
		-I/local/lib  \
		-I/usr/X11R6/include

VECOBJS= martin.o matpack.o sortdbl.o sort2dbl.o sort3dbl.o ivector.o \
		imatrix.o dvector.o dmatrix.o dmatmul.o cvector.o \
		cmatrix.o cmatmul.o dludcomp.o dsvdcomp.o dortho.o \
		cludcomp.o tred.o imtql.o imtql2.o rseigen.o cheigen.o \
		chgeigen.o matlapl.o moment.o matexph.o matsqth.o magicsquare.o

all : quartett.so

quartett.so: quartett.o   $(VECOBJS)
	$(CCC)  -Wall -o quartett.so $^  -shared -L$(MCPHASE_DIR)/bin/ic1ion_module/libs $(LFLAGS)
#-larpack -llapack -lblas -L/usr/local/gfortran/lib -lgfortran
# -static-libgcc -static-libgfortran  -L$(MCPHASE_DIR)/bin/ic1ion_module/libs -llapack -lblas -lgfortran

quartett.o: quartett.c
	$(CCC)  -Wall  -c  -fPIC -std=c++0x -o quartett.o quartett.c  $(INCLDIRS)

clean:
	rm -vf quartett.o 
	



