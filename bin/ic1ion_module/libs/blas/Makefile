#######################################################################
#  This is the minimum BLAS files required for ARPACK, from Netlib's
#  Reference BLAS implementation, 2007 version

include ../make.inc

#######################################################################
#  This is the makefile to create a library for the BLAS.
#  The files are grouped as follows:
#
#       SBLAS1 -- Single precision real BLAS routines
#       CBLAS1 -- Single precision complex BLAS routines
#       DBLAS1 -- Double precision real BLAS routines
#       ZBLAS1 -- Double precision complex BLAS routines
#
#       CB1AUX -- Real BLAS routines called by complex routines
#       ZB1AUX -- D.P. real BLAS routines called by d.p. complex
#                 routines
#
#      ALLBLAS -- Auxiliary routines for Level 2 and 3 BLAS
#
#       SBLAS2 -- Single precision real BLAS2 routines
#       CBLAS2 -- Single precision complex BLAS2 routines
#       DBLAS2 -- Double precision real BLAS2 routines
#       ZBLAS2 -- Double precision complex BLAS2 routines
#
#       SBLAS3 -- Single precision real BLAS3 routines
#       CBLAS3 -- Single precision complex BLAS3 routines
#       DBLAS3 -- Double precision real BLAS3 routines
#       ZBLAS3 -- Double precision complex BLAS3 routines
#
#  The library can be set up to include routines for any combination
#  of the four precisions.  To create or add to the library, enter make
#  followed by one or more of the precisions desired.  Some examples:
#       make single
#       make single complex
#       make single double complex complex16
#  Note that these commands are not safe for parallel builds.
#
#  Alternatively, the commands
#       make all
#  or
#       make
#  without any arguments creates a library of all four precisions.
#  The name of the library is held in BLASLIB, which is set in the
#  make.inc
#
#  To remove the object files after the library is created, enter
#       make clean
#  To force the source files to be recompiled, enter, for example,
#       make single FRC=FRC
#
#---------------------------------------------------------------------
#
#  Edward Anderson, University of Tennessee
#  March 26, 1990
#  Susan Ostrouchov,  September 30, 1994
#  Julie Langou, March 2007
#
#######################################################################

############################################################################
#  To create or add to the library, enter make followed by one or
#  more of the precisions desired.  Targets sdrv, ddrv, cdrv,
#  zdrv are used to add to the ARPACK library those BLAS routines needed by 
#  driver programs in the EXAMPLES directory.
#
#  Some examples:
#       make single
#       make single complex
#       make single sdrv
#       make single double complex complex16
#  Alternatively, the command
#       make
#  without any arguments creates a library of all four precisions.
#  The library also contains all extra BLAS routines used by driver
#  programs in the EXAMPLES directory.
#
#  The name of the library is defined by $(ARPACKLIB) in
#  ../ARmake.inc and is created at the next higher directory level.
#
#

SOBJ  = isamax.o sasum.o saxpy.o scopy.o sdot.o  snrm2.o srot.o  srotg.o\
        sscal.o  sswap.o sgemv.o sger.o  ssymv.o ssyr2.o strmm.o sgbmv.o

DOBJ  = idamax.o dasum.o daxpy.o dcopy.o ddot.o  dnrm2.o drot.o  drotg.o\
        dscal.o  dswap.o dgemv.o dger.o  dsymv.o dsyr2.o dtrmm.o dgbmv.o\
	dsyr2k.o dtrmv.o 

COBJ  = icamax.o caxpy.o ccopy.o cdotc.o cdotu.o cgbmv.o  cgemv.o cgerc.o\
        cgeru.o  cscal.o csscal.o cswap.o ctrmm.o ctrsv.o scasum.o scnrm2.o\
        isamax.o sscal.o 

ZOBJ  = izamax.o zaxpy.o zcopy.o zdotc.o zdotu.o zgbmv.o  zgemv.o zgerc.o\
        zgeru.o zscal.o zdscal.o zswap.o ztrmm.o ztrsv.o dzasum.o dznrm2.o\
        idamax.o dscal.o zhemm.o\
        zher2k.o zhemv.o zher2.o ztrmv.o dcabs1.o

# The following BLAS routine are not required by ARPACK source code.
# They are called by the driver programs in the EXAMPLES directory.
#
SDRVOBJ = sgemm.o strsm.o stbsv.o 
DDRVOBJ = dsymm.o dgemm.o dtrsm.o dtbsv.o
CDRVOBJ = ccopy.o cgemv.o cgeru.o cgemm.o cscal.o cswap.o ctrsm.o ctbsv.o\
          icamax.o 
ZDRVOBJ = zcopy.o zgemv.o zgeru.o zgemm.o zscal.o zswap.o ztrsm.o ztbsv.o\
          izamax.o

#---------------------------------------------------------------------

ALLOBJ = $(SOBJ) $(DOBJ) $(COBJ) $(ZOBJ) $(SDRVOBJ) $(DDRVOBJ) $(CDRVOBJ)\
         $(ZDRVOBJ)

#---------------------------------------------------------------------

#all: single double complex complex16 sdrv ddrv cdrv zdrv
all: $(BLASLIB)
 
$(BLASLIB): $(ALLOBJ)
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

single: $(SOBJ)
	$(ARCH) $(ARCHFLAGS) $(BLASLIB) $(SOBJ)
	$(RANLIB) $(BLASLIB)

double: $(DOBJ)
	$(ARCH) $(ARCHFLAGS) $(BLASLIB) $(DOBJ)
	$(RANLIB) $(BLASLIB)

complex: $(COBJ)
	$(ARCH) $(ARCHFLAGS) $(BLASLIB) $(COBJ)
	$(RANLIB) $(BLASLIB)

complex16: $(ZOBJ)
	$(ARCH) $(ARCHFLAGS) $(BLASLIB) $(ZOBJ)
	$(RANLIB) $(BLASLIB)
#
# Add routines needed by driver programs (in the EXAMPLES
# directory) to $(ARPACKLIB).
#
sdrv: $(SDRVOBJ) $(CDRVOBJ)
	$(ARCH) $(ARCHFLAGS) $(BLASLIB) $(SDRVOBJ) $(CDRVOBJ)
	$(RANLIB) $(BLASLIB)

ddrv: $(DDRVOBJ) $(ZDRVOBJ)
	$(ARCH) $(ARCHFLAGS) $(BLASLIB) $(DDRVOBJ) $(ZDRVOBJ)
	$(RANLIB) $(BLASLIB)

cdrv: $(CDRVOBJ)
	$(ARCH) $(ARCHFLAGS) $(BLASLIB) $(CDRVOBJ)
	$(RANLIB) $(BLASLIB)

zdrv: $(ZDRVOBJ)
	$(ARCH) $(ARCHFLAGS) $(BLASLIB) $(ZDRVOBJ)
	$(RANLIB) $(BLASLIB)

#---------------------------------------------------------------------

clean:
	rm -f *.o *.a

.f.o: 
	$(FORTRAN) $(OPTS) -c $< -o $@
