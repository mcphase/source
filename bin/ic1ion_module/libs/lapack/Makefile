############################################################################
#
#  Program:         ARPACK
#
#  Module:          Makefile
#
#  Purpose:         Sources Makefile
#
#  Creation date:   February 22, 1996
#
#  Modified:a       September 9, 1996
#
#  Send bug reports, comments or suggestions to arpack.caam.rice.edu
#
############################################################################
#\SCCS Information: @(#) 
# FILE: Makefile   SID: 2.2   DATE OF SID: 9/24/96   RELEASE: 2
 
include ../make.inc
 
############################################################################
#  To create or add to the library, enter make followed by one or
#  more of the precisions desired.  Targets sdrv, ddrv, cdrv,
#  zdrv are used to add to the library LAPACK routines needed by driver
#  programs in the EXAMPLES directory.
#
#  Some examples:
#       make single
#       make single sdrv
#       make single complex
#       make single double complex complex16
#  Alternatively, the command
#       make
#  without any arguments creates a library of all four precisions.
#  The library also contains all extra BLAS routines used by driver
#  programs in the EXAMPLES directory.
#
#  The name of the library is defined by $(LAPACKLIB) in
#  ../ARmake.inc and is created at the next higher directory level.
#
#
SOBJ  = sgeqr2.o slabad.o slacon.o slacpy.o sladiv.o slae2.o slaev2.o\
        slaexc.o slagtm.o slahqr.o slaln2.o slange.o slanhs.o\
        slanst.o slanv2.o slaptm.o slapy2.o slapy3.o slaran.o slarf.o\
        slarfg.o slarfx.o slarnd.o slarnv.o slartg.o slaruv.o slascl.o\
        slaset.o slasr.o  slasrt.o slassq.o slasy2.o sorm2r.o ssteqr.o\
        strevc.o strexc.o strsen.o strsyl.o

DOBJ  = dgeqr2.o dlabad.o dlacon.o dlacpy.o dladiv.o dlae2.o dlaev2.o\
        dlaexc.o dlagtm.o dlahqr.o dlaln2.o dlange.o dlanhs.o\
        dlanst.o dlanv2.o dlaptm.o dlapy2.o dlapy3.o dlaran.o dlarf.o\
        dlarfg.o dlarfx.o dlarnd.o dlarnv.o dlartg.o dlaruv.o dlascl.o\
        dlaset.o dlasr.o  dlasrt.o dlassq.o dlasy2.o dorm2r.o dsteqr.o\
        dtrevc.o dtrexc.o dtrsen.o dtrsyl.o \
	dsyev.o  dlansy.o dsytrd.o dsterf.o dorgtr.o dlatrd.o dsytd2.o\
	dorgqr.o dorgql.o dorg2r.o dlarft.o dlarfb.o dorg2l.o\
	dgesvd.o dgeqrf.o dorgbr.o dorglq.o dorgl2.o dgebrd.o dlabrd.o\
	dgebd2.o dbdsqr.o dlasv2.o dlas2.o  dlasq1.o dlasq2.o dlazq3.o\
	dlazq4.o dlasq5.o dlasq6.o dgelqf.o dgelq2.o dormbr.o dormqr.o\
	dormlq.o dorml2.o dlacn2.o \
	dlaebz.o dlagtf.o dlagts.o dlaneg.o dlarra.o dlarrb.o dlarrc.o\
	dlarrd.o dlarre.o dlarrf.o dlarrj.o dlarrk.o dlarrr.o dstebz.o\
	dlar1v.o dlarrv.o dorm2l.o dormql.o dormtr.o dstein.o dstemr.o\
	dstegr.o

IOBJ  = ilaenv.o lsame.o lsamen.o xerbla.o xlaenv.o iparmq.o ieeeck.o

CIOBJ = icmax1.o 

ZIOBJ = izmax1.o

COBJ  = cgeqr2.o clacon.o clacpy.o cladiv.o clahqr.o clange.o clanhs.o\
        clarf.o  clarfg.o clarnv.o clartg.o clascl.o claset.o classq.o\
        clatrs.o cmach.o  crot.o   ctrevc.o ctrexc.o ctrsen.o ctrsyl.o\
        cunm2r.o\
        scsum1.o slabad.o sladiv.o slamch.o slapy2.o slapy3.o slaruv.o

ZOBJ  = zgeqr2.o zlacon.o zlacpy.o zladiv.o zlahqr.o zlange.o zlanhs.o\
        zlarf.o  zlarfg.o zlarnv.o zlartg.o zlascl.o zlaset.o zlassq.o\
        zlatrs.o zmach.o  zrot.o   ztrevc.o ztrexc.o ztrsen.o ztrsyl.o\
        zunm2r.o\
        dzsum1.o dlabad.o dladiv.o dlamch.o dlapy2.o dlapy3.o dlaruv.o\
	zheev.o  zhesv.o  zhetf2.o zhetrf.o zhetrs.o zlahef.o zlanhe.o\
	zhetrd.o zungtr.o zsteqr.o zlatrd.o zhetd2.o zungql.o zungqr.o\
	zlasr.o  zung2l.o zlarft.o zlarfb.o zung2r.o\
	disnan.o dlaisnan.o\
	zlansy.o zlar1v.o zlarrv.o zstein.o zstemr.o zunm2l.o zunmql.o\
	zunmqr.o zunmtr.o zlacn2.o
#
# The following LAPACK routines are not needed by ARPACK
# source code.  They are used by the drivers in EXAMPLES
# directory.
#
SDRVOBJ = sgbtf2.o sgbtrf.o sgbtrs.o sgttrf.o sgttrs.o spttrf.o spttrs.o\
          slaswp.o sgtts2.o 

DDRVOBJ = dgbtf2.o dgbtrf.o dgbtrs.o dgttrf.o dgttrs.o dpttrf.o dpttrs.o\
          dlaswp.o dgtts2.o dsyevr.o

CDRVOBJ = cgbtf2.o cgbtrf.o cgbtrs.o cgttrf.o cgttrs.o claswp.o clacgv.o

ZDRVOBJ = zgbtrf.o zgbtrs.o zgbtf2.o zgttrf.o zgttrs.o zlaswp.o zlacgv.o\
          zheevr.o

#-----------------------------------------------------------------------
ALLOBJ=$(SOBJ) $(DOBJ) $(COBJ) $(ZOBJ) $(IOBJ) $(CIOBJ) $(ZIOBJ) \
       $(SDRVOBJ) $(DDRVOBJ) $(CDRVOBJ) $(ZDRVOBJ)
#-----------------------------------------------------------------------
GOTOOBJ = dsyev.o  dsteqr.o\
          dgeqr2.o\
	  dorm2r.o\
	  dlacpy.o dlarnv.o dlanst.o\
          dlapy2.o dlasr.o  dlartg.o dlascl.o dlaset.o dlaev2.o\
	  ilaenv.o\
	  dsytrd.o dsytd2.o dsterf.o\
	  dorgtr.o dorgqr.o dorgql.o dorg2r.o dorg2l.o\
	  dlae2.o  dlasrt.o dlarfg.o dlarf.o  dlaruv.o dlassq.o dlatrd.o dlarft.o dlarfb.o dlansy.o\
	  ieeeck.o iparmq.o
#-----------------------------------------------------------------------

.SUFFIXES:      .o .F .f

.f.o:
	$(FORTRAN) $(OPTS) -c $<

#all: single double complex complex16 sdrv ddrv cdrv zdrv
all: $(LAPACKLIB) libgotolapack.a

slamch.o: slamch.f ; $(FORTRAN) $(NOOPT) -c $< -o $@
dlamch.o: dlamch.f ; $(FORTRAN) $(NOOPT) -c $< -o $@

$(LAPACKLIB): $(ALLOBJ) slamch.o dlamch.o
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ) slamch.o dlamch.o
	$(RANLIB) $@

libgotolapack.a: $(GOTOOBJ) dlamch.o
	$(ARCH) $(ARCHFLAGS) $@ $(GOTOOBJ) dlamch.o
	$(RANLIB) $@

single: $(SOBJ) $(IOBJ)
	$(ARCH) $(ARCHFLAGS) $(LAPACKLIB) $(SOBJ) $(IOBJ)
	$(RANLIB) $(LAPACKLIB)
 
double: $(DOBJ) $(IOBJ)
	$(ARCH) $(ARCHFLAGS) $(LAPACKLIB) $(DOBJ) $(IOBJ)
	$(RANLIB) $(LAPACKLIB)

complex: $(COBJ) $(CIOBJ) $(IOBJ)
	$(ARCH) $(ARCHFLAGS) $(LAPACKLIB) $(COBJ) $(CIOBJ) $(IOBJ)
	$(RANLIB) $(LAPACKLIB)

complex16: $(ZOBJ) $(ZIOBJ) $(IOBJ)
	$(ARCH) $(ARCHFLAGS) $(LAPACKLIB) $(ZOBJ) $(ZIOBJ) $(IOBJ)
	$(RANLIB) $(LAPACKLIB) 
#
# Add routines needed by driver programs (in the EXAMPLES
# directory) to $(LAPACKLIB).
#
sdrv: $(SDRVOBJ) $(CDRVOBJ)
	$(ARCH) $(ARCHFLAGS) $(LAPACKLIB) $(SDRVOBJ) $(CDRVOBJ)
	$(RANLIB) $(LAPACKLIB)

ddrv: $(DDRVOBJ) $(ZDRVOBJ)
	$(ARCH) $(ARCHFLAGS) $(LAPACKLIB) $(DDRVOBJ) $(ZDRVOBJ)
	$(RANLIB) $(LAPACKLIB)

cdrv: $(CDRVOBJ)
	$(ARCH) $(ARCHFLAGS) $(LAPACKLIB) $(CDRVOBJ)
	$(RANLIB) $(LAPACKLIB)

zdrv: $(ZDRVOBJ)
	$(ARCH) $(ARCHFLAGS) $(LAPACKLIB) $(ZDRVOBJ)
	$(RANLIB) $(LAPACKLIB)

#  clean	- remove all object files
#
clean:
	rm -f *.o *.a
