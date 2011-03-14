# Makefile for all McPhase

include bin/src/Makefile.common

cfdir = bin/cf1ion_module
icdir = bin/ic1ion_module
mcpdir = bin/src
vecdir = bin/src/vector
funcdir = bin/src/functions
calcdir = bin/src/calc

all: vector functions cfield mcphase

vector: 
	cd $(vecdir) && $(MAKE) 

functions:
	cd $(funcdir) && $(MAKE)

calc:
	cd $(calcdir) && $(MAKE)

cfield: vector
	cd $(cfdir) && $(MAKE)

ic1ion: vector
	cd $(icdir) && $(MAKE)

mcphase: vector cfield ic1ion
	cd $(mcpdir) && $(MAKE)

clean: 
	cd $(vecdir) && $(MAKE) clean
	cd $(funcdir) && $(MAKE) clean
	cd $(cfdir) && $(MAKE) clean
	cd $(icdir) && $(MAKE) cleanall
	cd $(mcpdir) && $(MAKE) clean
	rm -vf bin/addj.exe bin/charges.exe bin/coq2jjj.exe \
		bin/mcdispit.exe bin/singleion.exe bin/cfield.exe \
		bin/cond.exe bin/jjj2j.exe bin/mcphasit.exe bin/spins.exe \
                bin/chrgplt.exe bin/pointc.exe bin/spinsfromq.exe \
                bin/mcdiff.exe bin/cf1ion_module/cfield.dll \
                bin/ic1ion.exe bin/icf1ion.exe bin/so1ion.exe \
                bin/ic1ion_module/ic1ion.dll \
                bin/fediff.exe bin/mf2fe.exe \
                bin/spindensplt.exe bin/orbmomdensplt.exe bin/momdensplt.exe \
                bin/spindensities.exe bin/orbmomdensities.exe bin/momdensities.exe \
                bin/currdensities.exe bin/currdensplt.exe \
                bin/formfactor.exe bin/radwavfunc.exe
	rm -vf bin/addj bin/charges bin/coq2jjj \
		bin/mcdispit bin/singleion bin/cfield \
		bin/cond bin/jjj2j bin/mcphasit bin/spins \
                bin/chrgplt bin/pointc bin/spinsfromq \
                bin/mcdiff bin/cf1ion_module/cfield.so \
                bin/ic1ion bin/icf1ion bin/so1ion \
                bin/ic1ion_module/ic1ion.so bin/ic1ion_module/icf1ion.so \
                bin/fediff bin/mf2fe \
                bin/spindensplt bin/orbmomdensplt bin/momdensplt \
                bin/spindensities bin/orbmomdensities bin/momdensities \
                bin/currdensities bin/currdensplt \
                bin/formfactor
