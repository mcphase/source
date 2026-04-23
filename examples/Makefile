#*********************************************************************
#
# File:          Makefile 
#
# Project:       McPhase
#
# Description:   Makefile for the examples
#
# Author(s):     M. Rotter
#
# Original:	 August 2001
#
#**********************************************************************

ifdef cross64
C64="cross64=1"
endif


vpath %.o  ../bin/src/vector

VECOBJS= matpack.o sortdbl.o sort2dbl.o sort3dbl.o ivector.o \
		imatrix.o dvector.o dmatrix.o dmatmul.o cvector.o \
		cmatrix.o cmatmul.o dludcomp.o dsvdcomp.o dortho.o \
		cludcomp.o tred.o imtql.o imtql2.o rseigen.o cheigen.o \
		chgeigen.o matlapl.o moment.o matexph.o matsqth.o magicsquare.o

all : cecu2a ndcu2b_1dim dyni2b2c erni2b2c gdni2b2c

FORCE:

cecu2a : FORCE $(VECOBJS) ./cecu2a/1ion_mod/kramer.c 
	cd ./cecu2a/1ion_mod ; make $(C64)

ndcu2b_1dim : FORCE   $(VECOBJS) ./ndcu2b_1dim/1ion_mod/kramer.c
	cd ./ndcu2b_1dim/1ion_mod ; make  $(C64)

dyni2b2c : FORCE  $(VECOBJS) ./dyni2b2c/1ion_mod/quartett.c
	cd ./dyni2b2c/1ion_mod ; make $(C64)

erni2b2c : FORCE  $(VECOBJS) ./erni2b2c/1ion_mod/quartett.c
	cd ./erni2b2c/1ion_mod ; make $(C64)

gdni2b2c : FORCE  $(VECOBJS) ./gdni2b2c/1ion_mod/brillouin.c  ./gdni2b2c/1ion_mod/bosefermi.c
	cd ./gdni2b2c/1ion_mod ; make $(C64)

clean :
	cd ./cecu2a/1ion_mod ; make clean
	cd ./ndcu2b_1dim/1ion_mod ; make clean
	cd ./dyni2b2c/1ion_mod ; make clean
	cd ./erni2b2c/1ion_mod ; make clean
	cd ./gdni2b2c/1ion_mod ; make clean
	rm -vf ./bfk/results/* ./bfk/results/.[!.]*
	cp ./.spins.eps ./bfk/results/
	cp ./readme.txt ./bfk/results/
	rm -vf ./cecu2a/results/* ./cecu2a/results/.[!.]*
	cp ./.spins.eps ./cecu2a/results/
	cp ./readme.txt ./cecu2a/results/
	rm -vf ./Ce3p_chain_cfphonon/results/* ./Ce3p_chain_cfphonon/results/.[!.]*
	cp ./.spins.eps ./Ce3p_chain_cfphonon/results/
	cp ./readme.txt ./Ce3p_chain_cfphonon/results/
	rm -vf ./Ce3p_tetragonalprim_cfphonon/results/* ./Ce3p_tetragonalprim_cfphonon/results/.[!.]*
	cp ./.spins.eps ./Ce3p_tetragonalprim_cfphonon/results/
	cp ./readme.txt ./Ce3p_tetragonalprim_cfphonon/results/
	rm -vf ./CeAl2_cfphonon_cfstrict/DMD_method/results/* ./CeAl2_cfphonon_cfstrict/CF/results/.[!.]*
	cp ./.spins.eps ./CeAl2_cfphonon_cfstrict/CF/results/
	cp ./readme.txt ./CeAl2_cfphonon_cfstrict/CF/results/
	rm -vf ./CeAl2_cfphonon_cfstrict/phonons/results/* ./CeAl2_cfphonon_cfstrict/phonons/results/.[!.]*
	cp ./.spins.eps ./CeAl2_cfphonon_cfstrict/phonons/results/
	cp ./readme.txt ./CeAl2_cfphonon_cfstrict/phonons/results/
	rm -vf ./cepd2si2/results/* ./cepd2si2/results/.[!.]*
	cp ./.spins.eps ./cepd2si2/results/
	cp ./readme.txt ./cepd2si2/results/
	rm -vf ./coo/results/* ./coo/results/.[!.]*
	cp ./.spins.eps ./coo/results/
	cp ./readme.txt ./coo/results/
	rm -vf ./dycu2iwata/results/* ./dycu2iwata/results/.[!.]*
	cp ./.spins.eps ./dycu2iwata/results/
	cp ./readme.txt ./dycu2iwata/results/
	rm -vf ./dyni2b2c/results/* ./dyni2b2c/results/.[!.]*
	cp ./.spins.eps ./dyni2b2c/results/
	cp ./readme.txt ./dyni2b2c/results/
	rm -vf ./erni2b2c/results/* ./erni2b2c/results/.[!.]*
	cp ./.spins.eps ./erni2b2c/results/
	cp ./readme.txt ./erni2b2c/results/
	rm -vf ./gd3gao6/results/* ./gd3gao6/results/.[!.]*
	cp ./.spins.eps ./gd3gao6/results/
	cp ./readme.txt ./gd3gao6/results/
	rm -vf ./gdni2b2c/results/* ./gdni2b2c/results/.[!.]*
	cp ./.spins.eps ./gdni2b2c/results/
	cp ./readme.txt ./gdni2b2c/results/
	rm -vf ./gdru2si2/results/* ./gdru2si2/results/.[!.]*
	cp ./.spins.eps ./gdru2si2/results/
	cp ./readme.txt ./gdru2si2/results/
	rm -vf ./helix_spinwave/results/* ./helix_spinwave/results/.[!.]*
	cp ./.spins.eps ./helix_spinwave/results/
	cp ./readme.txt ./helix_spinwave/results/
	rm -vf ./ho2ti2o7/results/* ./ho2ti2o7/results/.[!.]*
	cp ./.spins.eps ./ho2ti2o7/results/
	cp ./readme.txt ./ho2ti2o7/results/
	rm -vf ./la2coo4/results/* ./la2coo4/results/.[!.]*
	cp ./.spins.eps ./la2coo4/results/
	cp ./readme.txt ./la2coo4/results/
	rm -vf ./Cu2OSeO3/4interacting_tetrahedrons_noDM/results/* ./Cu2OSeO3/4interacting_tetrahedrons_noDM/results/.[!.]*
	cp ./readme.txt ./Cu2OSeO3/4interacting_tetrahedrons_noDM/results/
	rm -vf ./Cu2OSeO3/nocluster_spinwave_model_noDM/results/* ./Cu2OSeO3/nocluster_spinwave_model_noDM/results/.[!.]*
	cp ./readme.txt ./Cu2OSeO3/nocluster_spinwave_model_noDM/results/
	rm -vf ./LaCoO3_podlesnyak_polaron/2ions/results/* ./LaCoO3_podlesnyak_polaron/2ions/results/.[!.]*
	cp ./.spins.eps ./LaCoO3_podlesnyak_polaron/2ions/results/
	cp ./readme.txt ./LaCoO3_podlesnyak_polaron/2ions/results/
	rm -vf ./LaCoO3_podlesnyak_polaron/3ions/results/* ./LaCoO3_podlesnyak_polaron/3ions/results/.[!.]*
	cp ./.spins.eps ./LaCoO3_podlesnyak_polaron/3ions/results/
	cp ./readme.txt ./LaCoO3_podlesnyak_polaron/3ions/results/
	rm -vf ./LaCoO3_podlesnyak_polaron/5ions/results/* ./LaCoO3_podlesnyak_polaron/5ions/results/.[!.]*
	cp ./.spins.eps ./LaCoO3_podlesnyak_polaron/5ions/results/
	cp ./readme.txt ./LaCoO3_podlesnyak_polaron/5ions/results/
	rm -vf ./LaCoO3_podlesnyak_polaron/7ions/results/* ./LaCoO3_podlesnyak_polaron/7ions/results/.[!.]*
	cp ./.spins.eps ./LaCoO3_podlesnyak_polaron/7ions/results/
	cp ./readme.txt ./LaCoO3_podlesnyak_polaron/7ions/results/
	rm -vf ./lumno3/results/* ./lumno3/results/.[!.]*
	cp ./.spins.eps ./lumno3/results/
	cp ./readme.txt ./lumno3/results/
	rm -vf ./ndba2cu3o7/results/* ./ndba2cu3o7/results/.[!.]*
	cp ./.spins.eps ./ndba2cu3o7/results/
	cp ./readme.txt ./ndba2cu3o7/results/
	rm -vf ./ndcu2a/results/* ./ndcu2a/results/.[!.]*
	cp ./.spins.eps ./ndcu2a/results/
	cp ./readme.txt ./ndcu2a/results/
	rm -vf ./ndcu2b/results/* ./ndcu2b/results/.[!.]*
	cp ./.spins.eps ./ndcu2b/results/
	cp ./readme.txt ./ndcu2b/results/
	rm -vf ./ndcu2b_1dim/results/* ./ndcu2b_1dim/results/.[!.]*
	cp ./.spins.eps ./ndcu2b_1dim/results/
	cp ./readme.txt ./ndcu2b_1dim/results/
	rm -vf ./ndcu2b_new/results/* ./ndcu2b_new/results/.[!.]*
	cp ./.spins.eps ./ndcu2b_new/results/
	cp ./readme.txt ./ndcu2b_new/results/
	rm -vf ./ndcu2c/results/* ./ndcu2c/results/.[!.]*
	cp ./.spins.eps ./ndcu2c/results/
	cp ./readme.txt ./ndcu2c/results/
	rm -vf ./NiO/results/* ./NiO/results/.[!.]*
	cp ./.spins.eps ./NiO/results/
	cp ./readme.txt ./NiO/results/
	rm -vf ./Pr3Pd20Si6/results/* ./Pr3Pd20Si6/results/.[!.]*
	cp ./.spins.eps ./Pr3Pd20Si6/results/
	cp ./readme.txt ./Pr3Pd20Si6/results/
	rm -vf ./prni2b2c/results/* ./prni2b2c/results/.[!.]*
	cp ./.spins.eps ./prni2b2c/results/
	cp ./readme.txt ./prni2b2c/results/
	rm -rvf ./pupd3/results/* ./pupd3/results/.[!.]*
	cp ./.spins.eps ./pupd3/results/
	cp ./readme.txt ./pupd3/results/
	rm -vf ./Ru3p_create_sipf/results/* ./Ru3p_create_sipf/results/.[!.]*
	cp ./.spins.eps ./Ru3p_create_sipf/results/
	cp ./readme.txt ./Ru3p_create_sipf/results/
	rm -vf ./tbcu2a/results/* ./tbcu2a/results/.[!.]*
	cp ./.spins.eps ./tbcu2a/results/
	cp ./readme.txt ./tbcu2a/results/
	rm -vf ./prni2si2/results/* ./prni2si2/results/.[!.]*
	cp ./readme.txt ./prni2si2/results/
	rm -vf ./tmcu2_cf_phonon/results/* ./tmcu2_cf_phonon/results/.[!.]*
	rm -vf ./tmcu2_cf_phonon/mcphas.j
	cp ./readme.txt ./tmcu2_cf_phonon/results/
	rm -vf ./tmnic2/results/* ./tmnic2/results/.[!.]* ./tmnic2/mcphas.j ./tmnic2/mnn.j
	cp ./.spins.eps ./tmnic2/results/
	cp ./readme.txt ./tmnic2/results/
	rm -vf ./testic1ion/results/* ./testic1ion/results/.[!.]*
	cp ./readme.txt ./testic1ion/results/
	rm -vf ./tungsten_phonons/results/* ./tungsten_phonons/results/.[!.]*
	cp ./readme.txt ./tungsten_phonons/results/
	rm -vf ./upd3/results/mms/* ./upd3/results/mms/.[!.]*
	rm -vf ./upd3/results/* ./upd3/results/.[!.]*
	cp ./.spins.eps ./upd3/results/
	cp ./readme.txt ./upd3/results/
	rm -vf ../demo/results/* ../demo/results/.[!.]*
	cp ../examples/readme.txt ../demo/results/








