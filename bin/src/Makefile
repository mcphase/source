#*********************************************************************
#
# File:          Makefile
#
# Project:       McPhase
#
# Description:   Makefile for the program package - compiles and links all c++ code
#
# Author(s):     M. Rotter
#
# Original:	 August 2000
#
#**********************************************************************

# make looks in the following directories for files
vpath %.c  ./
vpath %.cpp ./
vpath %.h  ../include ../ic1ion_module
vpath %.hpp ../include ../ic1ion_module
vpath %.o  ./ ./vector ../cf1ion_module ./functions

#---------------------------------------------------------------------
#  Operating system
OSTYPE=$(shell uname)
ifdef cross64
OSTYPE=WIN
endif

#  (linux windows)
ifneq (,$(findstring $(OSTYPE),Linux))
 LINUX = linux
else
ifneq (,$(findstring $(OSTYPE),Darwin))
 LINUX = linux
else
 WINDOWS = windows
ifneq (,$(findstring CYGWIN,$(OSTYPE)))
 CYGWIN=1
endif
endif
endif
#
 STATIC = static
#
#---------------------------------------------------------------------
#               ... and his Flags
#
#                Includes
#

INCLDIRS =      -I. \
		-I../include \
		-I../ic1ion_module \
#		-I/usr/include/readline \
		-I/local/lib  \
		-I/usr/X11R6/include

ifdef CYGWIN
CYGFLAGS = -DCYGWIN
endif 

include ./Makefile.common

CFLAGS  = -c $(CCOPTS) $(INCLDIRS) $(CYGFLAGS) -U_FORTIFY_SOURCE

# PGPLLIBS =       -L$(PGPLOT_DIR)  -lcpgplot -lpgplot
ifdef LINUX
XLIBS= -ldl
#XLIBS=           -L/usr/X11/lib -lX11 -lXext -lc -lXpm -ldl
#EZLIB=           -ljpeg -ldl
else
EXE=.exe
endif
ifdef OMP
LFLAGS +=          $(XLIBS) $(EZLIB) -lstdc++ -lm -lgomp -lpthread
else
#LFLAGS =          $(XLIBS) $(EZLIB) -lstdc++ -lm
LFLAGS +=          $(XLIBS) $(EZLIB) -lstdc++ -lm ${LBLAS}
endif

# uncomment next line and set NUM_THREADS in mcphas_htcalc.c and mcdisp.c
# according to your machine if you want to use multiple processing to speed up
# calculations
 threads=1
# 
ifdef threads
THREADFLAG=-D_THREADS
ifdef LINUX
THREADLINK=-lpthread
endif
endif

#---------------------------------------------------------------------
# RCS options
#
#               RCS lock options
RCSLOCK  =      co -l -r$(VERSION)
#               RCS check out options
RCSCO    =      co -r$(VERSION)
#               RCS check in options
RCSCI    =      ci -u -f -s"Rel" -r$(VERSION) -m"$(LOCKMSG)"


#---------------------------------------------------------------------
# All files needed for the library
#
#	File list for RCS control
#               all include files

#               source files


#--------------------------------------------------------------------
# What has to be made
#
#               Names of the include files in the 
#		installation directories
#

INCL = mcphas.h  testspin.hpp qvectors.hpp spincf.hpp par.hpp  jjjpar.hpp inipar.hpp \
       mfcf.hpp vector.h martin.h graphic_parameters.hpp
OBJS = jjjpar.o par.o physprop.o qvectors.o testspin.o spincf.o inipar.o \
        mfcf.o martin.o ionpars.o graphic_parameters.o perlparse.o

VOBS = ivector.o imatrix.o dvector.o dmatrix.o dmatmul.o cvector.o cmatrix.o cmatmul.o \
	dludcomp.o dsvdcomp.o dortho.o cludcomp.o tred.o imtql.o imtql2.o rseigen.o cheigen.o \
	chgeigen.o matlapl.o moment.o matexph.o matsqth.o magicsquare.o \
	matpack.o sortdbl.o sort2dbl.o sort3dbl.o

CFIELDOBJS= cfieldrout.o diahermx.o eingabe.o intensit.o komplex.o matrix.o stevens.o  \
	theta.o ortho.o mb11a.o minima.o spline.o va05a.o

FUNCTIONOBJS=hypot3d.o hypot3f.o initds.o derf.o derfc.o ddaws.o ddaws2.o cwofz.o cwofz2.o \
	dcbrt.o dgamlm.o d9lgmc.o dlnrel.o dlngam.o dgamma.o dfac.o dbsi0e.o dbesi0.o dbsi1e.o \
	dbesi1.o dbsk0e.o dbesk0.o dbsk1e.o dbesk1.o d9b0mp.o d9b1mp.o dbesj0.o dbesj1.o dbesjn.o \
	dbesy0.o dbesy1.o dbesyn.o d9aimp.o daie.o dai.o dbie.o dbi.o ei.o jacobip.o chebyshevt.o \
	laguerrel.o hermiteh.o legendrep.o legendrepa.o harmonicy.o powi.o powii.o clnrel.o cgamma.o \
	clngam.o c9lgmc.o threejm.o threejj.o clebsch.o sixj.o

#LINKER= $(CCC) $(LDFLAGS)
#OBJS := $(patsubst %.c,%.o,$(wildcard ./src/*.c)) \
#	$(patsubst %.cpp,%.o,$(wildcard ./src/*.cpp)) 
#---------------------------------------------------------------------
# Build mcphas
#

all :  martin.o ../clusterize$(EXE) ../cond$(EXE) ../coq2jjj$(EXE) ../jjj2j$(EXE) \
	../singleion$(EXE) ../mcdispit$(EXE) ../mcphasit$(EXE) ../spins$(EXE) \
	../spinsfromq$(EXE) ../addj$(EXE)  ../reduce_unitcell$(EXE)  \
	../densplt$(EXE)  ../pointc$(EXE) ../mcdiff$(EXE) \
        ../formfactor$(EXE) ../radwavfunc$(EXE) ../anisotropyit$(EXE)

test: ../ionpars$(EXE) ../jjjpar$(EXE)
testperlparse:  ../perlparse$(EXE)
# ****************************************************
# compile sources 


martin.o:  martin.c martin.h vector.h
	$(CCC) $(CFLAGS) $< 

perlparse.o:  perlparse.c martin.h vector.h myev.h perlparse.h
	$(CCC) $(CFLAGS) $< 

graphic_parameters.o: graphic_parameters.cpp graphic_parameters.hpp martin.h vector.h
	$(CCC) $(CFLAGS) $<

paramete.o: paramete.cpp paramete.hpp martin.h vector.h 
	$(CCC) $(CFLAGS) $< 

par.o: par.cpp par.hpp jjjpar.hpp martin.h vector.h ionpars.hpp cryststruct.hpp
	$(CCC) $(CFLAGS) $< 

jjjpar.o: jjjpar.cpp par.hpp ionpars.hpp perlparse.h jjjpar_observables.cpp jjjpar_intmod_kramer.cpp jjjpar_intmod_brillouin.cpp jjjpar_intmod_cluster.cpp jjjpar_basmodfunc.cpp jjjpar.hpp martin.h vector.h myev.h
	$(CCC) $(CFLAGS) $< 

testspin.o: testspin.cpp testspin.hpp spincf.hpp martin.h
	$(CCC) $(CFLAGS) $< 

spincf.o: spincf.cpp spincf_out.cpp spincf.hpp vector.h density.hpp jjjpar.hpp graphic_parameters.hpp cryststruct.hpp
	$(CCC) $(CFLAGS) $< 

density.o: density.cpp density.hpp vector.h
	$(CCC) $(CFLAGS) $< 


mfcf.o: mfcf.cpp mfcf.hpp vector.h
	$(CCC) $(CFLAGS) $< 

mdcf.o: mdcf.cpp mdcf.hpp vector.h
	$(CCC) $(CFLAGS) $< 

jq.o: jq.cpp jq.hpp vector.h mdcf.hpp
	$(CCC) $(CFLAGS) $< 

inipar.o: inipar.cpp inipar.hpp martin.h vector.h par.hpp
	$(CCC) $(CFLAGS) $< 

ionpars.o: ionpars.cpp ionpars.hpp perlparse.h martin.h
	$(CCC) $(CFLAGS) $<

inimcdis.o: inimcdis.cpp inimcdis.hpp martin.h vector.h mfcf.hpp 
	$(CCC) $(CFLAGS) $< 

physprop.o: physprop.cpp physprop.hpp spincf.hpp vector.h mfcf.hpp martin.h par.hpp inipar.hpp
	$(CCC) $(CFLAGS) $< 

qvectors.o: qvectors.cpp qvectors.hpp spincf.hpp martin.h vector.h inipar.hpp par.hpp myev.c
	$(CCC) $(CFLAGS) $<

mcphas.o: mcphas.c $(INCL) mcphas_htcalc.c mcphas_fecalc.c mcphas_physpropcalc.c
#	$(CCC) $(CFLAGS)  $<
	$(CCC) $(CFLAGS) $(THREADFLAG) $<

anisotropy.o: anisotropy.c $(INCL) mcphas_htcalc.c mcphas_fecalc.c mcphas_physpropcalc.c
#	$(CCC) $(CFLAGS)   $<
	$(CCC) $(CFLAGS)  $(THREADFLAG) $<


mcdisp.o: mcdisp.c  myev.c mcdisp_intcalc.c mcdisp_output.c jq.hpp mdcf.hpp trs_io.c $(INCL)
#	$(CCC) $(CFLAGS)  $<
	$(CCC) $(CFLAGS) $(THREADFLAG) $<


cond.o: cond.c paramete.hpp martin.h
	$(CCC) $(CFLAGS) $<

clusterize.o: clusterize.c par.hpp martin.h jjjpar.hpp 
	$(CCC) $(CFLAGS) $<

coq2jjj.o: coq2jjj.c paramete.hpp martin.h
	$(CCC) $(CFLAGS) $<

jjj2j.o: jjj2j.c par.hpp jjj2j.h jjjpar.hpp vector.h
	$(CCC) $(CFLAGS) $<

spins.o: spins.c densities_func.c spincf.hpp  martin.h vector.h graphic_parameters.hpp cryststruct.hpp par.hpp jjjpar.hpp
	$(CCC) $(CFLAGS) $<

spinsfromq.o: spinsfromq.c spincf.hpp  par.hpp martin.h vector.h 
	$(CCC) $(CFLAGS) $<

singleion.o: singleion.c par.hpp jjjpar.hpp myev.h martin.h vector.h trs_io.c
	$(CCC) $(CFLAGS) $<

densplt.o:  densplt.c plt_func.c density.hpp martin.h vector.h jjjpar.hpp graphic_parameters.hpp cryststruct.hpp spincf.hpp
	$(CCC) $(CFLAGS) $<

pointc.o:  pointc.c  martin.h vector.h ionpars.hpp jjjpar.hpp
	$(CCC) $(CFLAGS) $<

formfactor.o:  formfactor.c  martin.h vector.h ionpars.hpp jjjpar.hpp
	$(CCC) $(CFLAGS) $<

radwavfunc.o:  radwavfunc.c  martin.h vector.h ionpars.hpp jjjpar.hpp
	$(CCC) $(CFLAGS) $<

addj.o: addj.c par.hpp  martin.h vector.h
	$(CCC) $(CFLAGS) $<

reduce_unitcell.o: reduce_unitcell.c par.hpp  martin.h vector.h
	$(CCC) $(CFLAGS) $<

mcdiff.o: mcdiff.c jjjpar.hpp martin.h myev.h mcdiff.h mcdiff_intcalc.c mcdiff_output.c spincf.hpp mfcf.hpp
	$(CCC) $(CFLAGS) $<

myev.o: myev.c
	$(CCC) $(CFLAGS) $<

mf2fe.o: mf2fe.c
	$(CCC) $(CFLAGS) $<

fediff.o: fediff.c
	$(CCC) $(CFLAGS) $<

ifdef LINUX

inimenu.o: inimenu.c
	$(CCC) $(CFLAGS) $<

endif

ip.o: ionpars.cpp
	$(CCC) $(CFLAGS) -o  $<

# ************************************************
# make  binaries (link .o files)

#for testing only
../ionpars$(EXE): ionpars.o $(CFIELDOBJS) $(VOBS) $(FUNCTIONOBJS) martin.o myev.o
	@echo Linking ../ionpars
	@$(LINKER)  -Wall  -o ../ionpars  $^ $(LFLAGS)

../jjjpar$(EXE): jjjpar.o ionpars.o $(CFIELDOBJS) $(VOBS) $(FUNCTIONOBJS) martin.o myev.o
	@echo Linking ../jjjpar
	@$(LINKER)  -Wall  -o ../jjjpar  $^ $(LFLAGS)

../perlparse$(EXE): perlparse.o $(CFIELDOBJS) $(VOBS) $(FUNCTIONOBJS) martin.o myev.o
	@echo Linking ../perlparse
	@$(LINKER)  -Wall  -o ../perlparse  $^ $(LFLAGS)

# mcphase software suite programs
../cond$(EXE): cond.o paramete.o $(VOBS) martin.o
	@echo Linking ../cond
	@$(LINKER)  -Wall  -o ../cond $^ $(LFLAGS)

../clusterize$(EXE): clusterize.o $(VOBS) martin.o par.o jjjpar.o ionpars.o $(FUNCTIONOBJS) myev.o  $(CFIELDOBJS)   perlparse.o
	@echo Linking ../clusterize
	@$(LINKER)  -Wall  -o ../clusterize $^ $(LFLAGS)

../coq2jjj$(EXE): coq2jjj.o paramete.o $(VOBS) martin.o
	@echo Linking ../coq2jjj
	@$(LINKER)  -Wall  -o ../coq2jjj $^ $(LFLAGS)

../jjj2j$(EXE): jjj2j.o par.o jjjpar.o  $(VOBS) martin.o $(CFIELDOBJS) $(FUNCTIONOBJS) ionpars.o myev.o perlparse.o
	@echo Linking ../jjj2j
	@$(LINKER)  -Wall  -o ../jjj2j $^ $(LFLAGS)

../mcdispit$(EXE): mcdisp.o martin.o par.o jjjpar.o inimcdis.o jq.o mdcf.o mfcf.o ionpars.o perlparse.o $(VOBS) $(FUNCTIONOBJS)  $(CFIELDOBJS)
	@echo Linking ../mcdispit
	@$(LINKER)  -Wall  -o ../mcdispit   $^  $(LFLAGS) $(THREADLINK)

../mcphasit$(EXE): $(OBJS) martin.o $(VOBS)  $(CFIELDOBJS) $(FUNCTIONOBJS) ionpars.o density.o mcphas.o
	@echo Linking ../mcphasit
	@$(LINKER)  -Wall  -o ../mcphasit  $^ $(LFLAGS) $(THREADLINK)

../anisotropyit$(EXE): $(OBJS) martin.o $(VOBS)  $(CFIELDOBJS) $(FUNCTIONOBJS) ionpars.o perlparse.o density.o anisotropy.o
	@echo Linking ../anisotropyit
	@$(LINKER)  -Wall  -o ../anisotropyit  $^ $(LFLAGS) $(THREADLINK)

../spins$(EXE): spincf.o spins.o martin.o myev.o  $(VOBS) density.o ionpars.o jjjpar.o perlparse.o $(CFIELDOBJS) $(FUNCTIONOBJS) graphic_parameters.o par.o
	@echo Linking ../spins
	@$(LINKER)  -Wall  -o ../spins $^  $(LFLAGS)

../spinsfromq$(EXE): spincf.o spinsfromq.o martin.o jjjpar.o par.o myev.o perlparse.o $(FUNCTIONOBJS)  $(VOBS) $(CFIELDOBJS) ionpars.o density.o
	@echo Linking ../spinsfromq
	@$(LINKER)  -Wall  -o ../spinsfromq $^  $(LFLAGS)

../singleion$(EXE):  singleion.o martin.o  jjjpar.o par.o myev.o perlparse.o $(VOBS) $(FUNCTIONOBJS)  $(CFIELDOBJS) ionpars.o density.o
	@echo Linking ../singleion
	@$(LINKER)  -Wall  -o ../singleion $^  $(LFLAGS)

../densplt$(EXE): densplt.o martin.o $(VOBS) $(CFIELDOBJS) myev.o perlparse.o $(FUNCTIONOBJS) ionpars.o density.o spincf.o jjjpar.o graphic_parameters.o  par.o
	@echo Linking ../densplt
	@$(LINKER)  -Wall  -o ../densplt $^  $(LFLAGS)

../pointc$(EXE): pointc.o martin.o $(VOBS) $(CFIELDOBJS) ionpars.o myev.o perlparse.o  $(FUNCTIONOBJS) jjjpar.o  par.o
	@echo Linking ../pointc
	@$(LINKER)  -Wall  -o ../pointc $^  $(LFLAGS)

../formfactor$(EXE): formfactor.o martin.o $(VOBS) $(CFIELDOBJS) ionpars.o myev.o perlparse.o $(FUNCTIONOBJS) jjjpar.o  par.o
	@echo Linking ../formfactor
	@$(LINKER)  -Wall  -o ../formfactor $^  $(LFLAGS)

../radwavfunc$(EXE): radwavfunc.o martin.o $(VOBS) $(CFIELDOBJS) ionpars.o myev.o perlparse.o  $(FUNCTIONOBJS) jjjpar.o  par.o
	@echo Linking ../radwavfunc
	@$(LINKER)  -Wall  -o ../radwavfunc $^  $(LFLAGS)

../addj$(EXE): addj.o par.o martin.o myev.o jjjpar.o perlparse.o par.o $(VOBS) $(CFIELDOBJS) $(FUNCTIONOBJS)  ionpars.o density.o
	@echo Linking ../addj
	@$(LINKER)  -Wall  -o ../addj $^  $(LFLAGS)

../reduce_unitcell$(EXE): reduce_unitcell.o par.o martin.o myev.o jjjpar.o perlparse.o $(VOBS) $(CFIELDOBJS) $(FUNCTIONOBJS)  ionpars.o density.o
	@echo Linking ../reduce_unitcell
	@$(LINKER)  -Wall  -o ../reduce_unitcell $^  $(LFLAGS)


../mcdiff$(EXE): mcdiff.o jjjpar.o martin.o myev.o perlparse.o $(VOBS) $(FUNCTIONOBJS) $(CFIELDOBJS) ionpars.o density.o spincf.o mfcf.o par.o
	@echo Linking ../mcdiff
	@$(LINKER)  -Wall  -o ../mcdiff $^  $(LFLAGS)

../mf2fe$(EXE): mf2fe.o martin.o jjjpar.o inipar.o par.o ionpars.o density.o spincf.o mfcf.o testspin.o myev.o $(VOBS) $(CFIELDOBJS) $(FUNCTIONOBJS)
	@echo Linking ../mf2fe
	@$(LINKER)  -Wall  -o ../mf2fe $^  $(LFLAGS)

../fediff$(EXE): fediff.o martin.o jjjpar.o inipar.o par.o ionpars.o density.o spincf.o mfcf.o testspin.o myev.o $(VOBS) $(CFIELDOBJS) $(FUNCTIONOBJS)
	@echo Linking ../fediff
	@$(LINKER)  -Wall  -o ../fediff $^  $(LFLAGS)

ifdef LINUX

inimenu: inimenu.o inipar.o martin.o $(VOBS)
	@echo Linking ../inimenu
	@$(LINKER)  -Wall  -o ../inimenu $^  $(LFLAGS)

endif

clean: 
	rm -vf ./*.o

