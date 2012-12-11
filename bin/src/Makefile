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
vpath %.h  ../include
vpath %.hpp ../include
vpath %.o  ./ ./vector ../cf1ion_module ./functions

#---------------------------------------------------------------------
#  Operating system
OSTYPE=$(shell uname)
#  (linux windows)
ifneq (,$(findstring $(OSTYPE),Linux))
 LINUX = linux
else
ifneq (,$(findstring $(OSTYPE),Darwin))
 LINUX = linux
else
 WINDOWS = windows
ifneq (,$(findstring CYGWIN,$(shell uname)))
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
OBJS = jjjpar.o par.o physprop.o qvectors.o mcphas.o testspin.o spincf.o inipar.o \
        mfcf.o martin.o ionpars.o graphic_parameters.o

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

LINKER= $(CCC) $(LDFLAGS)
#OBJS := $(patsubst %.c,%.o,$(wildcard ./src/*.c)) \
#	$(patsubst %.cpp,%.o,$(wildcard ./src/*.cpp)) 
#---------------------------------------------------------------------
# Build mcphas
#

all :  martin.o ../cond$(EXE) ../coq2jjj$(EXE) ../jjj2j$(EXE) \
	../mcdispit$(EXE) ../mcphasit$(EXE) ../spins$(EXE) \
	../spinsfromq$(EXE) ../addj$(EXE) ../singleion$(EXE) ../charges$(EXE)  \
	../chrgplt$(EXE) ../spindensplt$(EXE) ../pointc$(EXE) ../mcdiff$(EXE) \
        ../orbmomdensplt$(EXE) ../momdensplt$(EXE) ../currdensplt$(EXE) \
        ../formfactor$(EXE) ../radwavfunc$(EXE) ../spindensities$(EXE) \
        ../orbmomdensities$(EXE) ../momdensities$(EXE) ../currdensities$(EXE)



test: ../ionpars$(EXE) ../jjjpar$(EXE)
# ****************************************************
# compile sources 


martin.o:  martin.c martin.h vector.h
	$(CCC) $(CFLAGS) $< 
	
graphic_parameters.o: graphic_parameters.cpp graphic_parameters.hpp martin.h vector.h
	$(CCC) $(CFLAGS) $<

paramete.o: paramete.cpp paramete.hpp martin.h vector.h 
	$(CCC) $(CFLAGS) $< 

par.o: par.cpp par.hpp jjjpar.hpp martin.h vector.h
	$(CCC) $(CFLAGS) $< 

jjjpar.o: jjjpar.cpp jjjpar_intmod_kramer.cpp jjjpar_intmod_brillouin.cpp jjjpar_intmod_cluster.cpp jjjpar_basmodfunc.cpp jjjpar.hpp martin.h vector.h myev.h
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

jq.o: jq.cpp jq.hpp vector.h
	$(CCC) $(CFLAGS) $< 

inipar.o: inipar.cpp inipar.hpp martin.h vector.h 
	$(CCC) $(CFLAGS) $< 

ionpars.o: ionpars.cpp ionpars.hpp
	$(CCC) $(CFLAGS) $<

inimcdis.o: inimcdis.cpp inimcdis.hpp martin.h vector.h mfcf.hpp
	$(CCC) $(CFLAGS) $< 

physprop.o: physprop.cpp physprop.hpp spincf.hpp vector.h mfcf.hpp martin.h par.hpp
	$(CCC) $(CFLAGS) $< 

qvectors.o: qvectors.cpp qvectors.hpp spincf.hpp martin.h vector.h inipar.hpp myev.c
	$(CCC) $(CFLAGS) $<

mcphas.o: mcphas.c $(INCL) mcphas_htcalc.c mcphas_fecalc.c mcphas_physpropcalc.c
	$(CCC) $(CFLAGS)  $<

mcdisp.o: mcdisp.c  myev.c mcdisp_intcalc.c jq.hpp mdcf.hpp $(INCL)
	$(CCC) $(CFLAGS) $(THREADFLAG) $<

cond.o: cond.c paramete.hpp martin.h
	$(CCC) $(CFLAGS) $<

coq2jjj.o: coq2jjj.c paramete.hpp martin.h
	$(CCC) $(CFLAGS) $<

jjj2j.o: jjj2j.c par.hpp jjj2j.h jjjpar.hpp vector.h
	$(CCC) $(CFLAGS) $<

spins.o: spins.c densities_func.c spincf.hpp  martin.h vector.h graphic_parameters.hpp cryststruct.hpp
	$(CCC) $(CFLAGS) $<

spinsfromq.o: spinsfromq.c spincf.hpp  par.hpp martin.h vector.h 
	$(CCC) $(CFLAGS) $<

singleion.o: singleion.c par.hpp  martin.h vector.h
	$(CCC) $(CFLAGS) $<

charges.o: charges.c densities_func.c par.hpp  spincf.hpp martin.h vector.h myev.h graphic_parameters.hpp cryststruct.hpp
	$(CCC) $(CFLAGS) $<

spindensities.o: spindensities.c densities_func.c par.hpp  spincf.hpp martin.h vector.h myev.h graphic_parameters.hpp cryststruct.hpp
	$(CCC) $(CFLAGS) $<

orbmomdensities.o: orbmomdensities.c densities_func.c par.hpp  spincf.hpp martin.h vector.h myev.h graphic_parameters.hpp cryststruct.hpp
	$(CCC) $(CFLAGS) $<

momdensities.o: momdensities.c densities_func.c par.hpp  spincf.hpp martin.h vector.h myev.h graphic_parameters.hpp cryststruct.hpp
	$(CCC) $(CFLAGS) $<

currdensities.o: currdensities.c densities_func.c par.hpp  spincf.hpp martin.h vector.h myev.h graphic_parameters.hpp cryststruct.hpp
	$(CCC) $(CFLAGS) $<

chrgplt.o:  chrgplt.c plt_func.c density.hpp martin.h vector.h jjjpar.hpp graphic_parameters.hpp cryststruct.hpp
	$(CCC) $(CFLAGS) $<

spindensplt.o:  spindensplt.c plt_func.c density.hpp martin.h vector.h jjjpar.hpp graphic_parameters.hpp cryststruct.hpp
	$(CCC) $(CFLAGS) $<

orbmomdensplt.o:  orbmomdensplt.c plt_func.c density.hpp martin.h vector.h jjjpar.hpp graphic_parameters.hpp cryststruct.hpp
	$(CCC) $(CFLAGS) $<

momdensplt.o:  momdensplt.c plt_func.c density.hpp martin.h vector.h jjjpar.hpp graphic_parameters.hpp cryststruct.hpp
	$(CCC) $(CFLAGS) $<

currdensplt.o:  currdensplt.c plt_func.c density.hpp martin.h vector.h jjjpar.hpp graphic_parameters.hpp cryststruct.hpp
	$(CCC) $(CFLAGS) $<


pointc.o:  pointc.c  martin.h vector.h ionpars.hpp jjjpar.hpp
	$(CCC) $(CFLAGS) $<

formfactor.o:  formfactor.c  martin.h vector.h ionpars.hpp jjjpar.hpp
	$(CCC) $(CFLAGS) $<

radwavfunc.o:  radwavfunc.c  martin.h vector.h ionpars.hpp jjjpar.hpp
	$(CCC) $(CFLAGS) $<

addj.o: addj.c par.hpp  martin.h vector.h
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

# mcphase software suite programs
../cond$(EXE): cond.o paramete.o $(VOBS) martin.o
	@echo Linking ../cond
	@$(LINKER)  -Wall  -o ../cond $^ $(LFLAGS)

../coq2jjj$(EXE): coq2jjj.o paramete.o $(VOBS) martin.o
	@echo Linking ../coq2jjj
	@$(LINKER)  -Wall  -o ../coq2jjj $^ $(LFLAGS)

../jjj2j$(EXE): jjj2j.o par.o jjjpar.o  $(VOBS) martin.o $(CFIELDOBJS) $(FUNCTIONOBJS) ionpars.o myev.o
	@echo Linking ../jjj2j
	@$(LINKER)  -Wall  -o ../jjj2j $^ $(LFLAGS)

../mcdispit$(EXE): mcdisp.o martin.o par.o jjjpar.o inimcdis.o jq.o mdcf.o mfcf.o ionpars.o  $(VOBS) $(FUNCTIONOBJS)  $(CFIELDOBJS)
	@echo Linking ../mcdispit
	@$(LINKER)  -Wall  -o ../mcdispit   $^  $(LFLAGS) $(THREADLINK)

../mcphasit$(EXE): $(OBJS) martin.o $(VOBS)  $(CFIELDOBJS) $(FUNCTIONOBJS) ionpars.o density.o
	@echo Linking ../mcphasit
	@$(LINKER)  -Wall  -o ../mcphasit  $^ $(LFLAGS) 

../spins$(EXE): spincf.o spins.o martin.o myev.o  $(VOBS) density.o ionpars.o jjjpar.o $(CFIELDOBJS) $(FUNCTIONOBJS) graphic_parameters.o par.o
	@echo Linking ../spins
	@$(LINKER)  -Wall  -o ../spins $^  $(LFLAGS)

../spinsfromq$(EXE): spincf.o spinsfromq.o martin.o jjjpar.o par.o myev.o $(FUNCTIONOBJS)  $(VOBS) $(CFIELDOBJS) ionpars.o density.o
	@echo Linking ../spinsfromq
	@$(LINKER)  -Wall  -o ../spinsfromq $^  $(LFLAGS)

../singleion$(EXE):  singleion.o martin.o  jjjpar.o par.o myev.o $(VOBS) $(FUNCTIONOBJS)  $(CFIELDOBJS) ionpars.o density.o
	@echo Linking ../singleion
	@$(LINKER)  -Wall  -o ../singleion $^  $(LFLAGS)

../charges$(EXE):  charges.o martin.o  jjjpar.o par.o $(VOBS) $(FUNCTIONOBJS) $(CFIELDOBJS) ionpars.o spincf.o density.o myev.o graphic_parameters.o  par.o
	@echo Linking ../charges
	@$(LINKER)  -Wall  -o ../charges $^  $(LFLAGS)

../spindensities$(EXE):  spindensities.o martin.o  jjjpar.o par.o $(VOBS) $(FUNCTIONOBJS) $(CFIELDOBJS) ionpars.o spincf.o density.o myev.o graphic_parameters.o  par.o
	@echo Linking ../spindensities
	@$(LINKER)  -Wall  -o ../spindensities $^  $(LFLAGS)

../orbmomdensities$(EXE):  orbmomdensities.o martin.o  jjjpar.o par.o $(VOBS) $(FUNCTIONOBJS) $(CFIELDOBJS) ionpars.o spincf.o density.o myev.o graphic_parameters.o  par.o
	@echo Linking ../orbmomdensities
	@$(LINKER)  -Wall  -o ../orbmomdensities $^  $(LFLAGS)

../momdensities$(EXE):  momdensities.o martin.o  jjjpar.o par.o $(VOBS) $(FUNCTIONOBJS) $(CFIELDOBJS) ionpars.o spincf.o density.o myev.o graphic_parameters.o  par.o
	@echo Linking ../momdensities
	@$(LINKER)  -Wall  -o ../momdensities $^  $(LFLAGS)

../currdensities$(EXE):  currdensities.o martin.o  jjjpar.o par.o $(VOBS) $(FUNCTIONOBJS) $(CFIELDOBJS) ionpars.o spincf.o density.o myev.o graphic_parameters.o  par.o
	@echo Linking ../currdensities
	@$(LINKER)  -Wall  -o ../currdensities $^  $(LFLAGS)

../chrgplt$(EXE): chrgplt.o martin.o $(VOBS) $(CFIELDOBJS) myev.o $(FUNCTIONOBJS) ionpars.o density.o spincf.o jjjpar.o graphic_parameters.o  par.o
	@echo Linking ../chrgplt
	@$(LINKER)  -Wall  -o ../chrgplt $^  $(LFLAGS)

../spindensplt$(EXE): spindensplt.o martin.o $(VOBS) $(CFIELDOBJS) myev.o  $(FUNCTIONOBJS) ionpars.o density.o spincf.o jjjpar.o graphic_parameters.o  par.o
	@echo Linking ../spindensplt
	@$(LINKER)  -Wall  -o ../spindensplt $^  $(LFLAGS)

../orbmomdensplt$(EXE): orbmomdensplt.o martin.o $(VOBS) myev.o  $(CFIELDOBJS) $(FUNCTIONOBJS) ionpars.o density.o spincf.o jjjpar.o graphic_parameters.o  par.o
	@echo Linking ../orbmomdensplt
	@$(LINKER)  -Wall  -o ../orbmomdensplt $^  $(LFLAGS)

../momdensplt$(EXE): momdensplt.o martin.o $(VOBS) $(CFIELDOBJS) myev.o  $(FUNCTIONOBJS) ionpars.o density.o spincf.o jjjpar.o graphic_parameters.o  par.o
	@echo Linking ../momdensplt
	@$(LINKER)  -Wall  -o ../momdensplt $^  $(LFLAGS)

../currdensplt$(EXE): currdensplt.o martin.o $(VOBS) $(CFIELDOBJS) myev.o  $(FUNCTIONOBJS) ionpars.o density.o spincf.o jjjpar.o graphic_parameters.o  par.o
	@echo Linking ../currdensplt
	@$(LINKER)  -Wall  -o ../currdensplt $^  $(LFLAGS)

../pointc$(EXE): pointc.o martin.o $(VOBS) $(CFIELDOBJS) ionpars.o myev.o  $(FUNCTIONOBJS) jjjpar.o  par.o
	@echo Linking ../pointc
	@$(LINKER)  -Wall  -o ../pointc $^  $(LFLAGS)

../formfactor$(EXE): formfactor.o martin.o $(VOBS) $(CFIELDOBJS) ionpars.o myev.o  $(FUNCTIONOBJS) jjjpar.o  par.o
	@echo Linking ../formfactor
	@$(LINKER)  -Wall  -o ../formfactor $^  $(LFLAGS)

../radwavfunc$(EXE): radwavfunc.o martin.o $(VOBS) $(CFIELDOBJS) ionpars.o myev.o  $(FUNCTIONOBJS) jjjpar.o  par.o
	@echo Linking ../radwavfunc
	@$(LINKER)  -Wall  -o ../radwavfunc $^  $(LFLAGS)

../addj$(EXE): addj.o par.o martin.o myev.o jjjpar.o $(VOBS) $(CFIELDOBJS) $(FUNCTIONOBJS)  ionpars.o density.o
	@echo Linking ../addj
	@$(LINKER)  -Wall  -o ../addj $^  $(LFLAGS)

../mcdiff$(EXE): mcdiff.o jjjpar.o martin.o myev.o $(VOBS) $(FUNCTIONOBJS) $(CFIELDOBJS) ionpars.o density.o spincf.o mfcf.o par.o
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

