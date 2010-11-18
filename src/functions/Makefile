#################################################################################
#  										#
# Makefile for the MatPack 1.0 special functions library  by B.M.Gammel, 1995	#
#  										#
# Tested with the following C++ compilers and architectures:			#
#  										#
#   GNU C++ 2.5.8 compiler of the Free Software Foundation			#
#       on  SUN Sparc     (SunOs 4.1.2)						#
#           SUN SLC 2/10  (SunOs 4.1.2)						#
#           DEC 3100/5000 (Ultrix 4.2)						#
#  										#
#   GNU C++ 2.6.3 on Linux Pentium PC						#
#	    2.7.2 on DEC Alpha OSF 1.3A 					#
#  										#
#   DEC C++ cxx 1.3A compiler	on  DEC Alpha OSF 1.3A				#
#  										#
# Permission to use , copy, and  distribute  Matpack  in  its entirety  and its #
# documentation  for non-commercial purpose and  without fee is hereby granted, #
# provided that this license information and copyright notice appear unmodified #
# in all copies.  This software is provided 'as is'  without express or implied #
# warranty.  In no event will the author be held liable for any damages arising #
# from the use of this software.						#
# Note that distributing Matpack 'bundled' in with any product is considered to #
# be a 'commercial purpose'.							#
# The software may be modified for your own purposes, but modified versions may #
# not be distributed without prior consent of the author.			#
#                                                                               #
# Read the  COPYRIGHT and  README files in this distribution about registration	#
# and installation of Matpack.							#
#                                                                               #
#################################################################################

###########################################################################
### set path to matpack home directory
###########################################################################

MATPACK = ..

###########################################################################
### include common rules
###########################################################################

include ${MATPACK}/Makefile.common


###########################################################################
### make all objects
###########################################################################

FCNINC= ../../include/mpspecfunp.h

FCNOBJ= hypot3d.o  hypot3f.o  \
	initds.o  derf.o  derfc.o  ddaws.o  ddaws2.o  cwofz.o  cwofz2.o  \
	dcbrt.o  \
	dgamlm.o  d9lgmc.o  dlnrel.o  \
	dlngam.o  dgamma.o  dfac.o  \
	dbsi0e.o  dbesi0.o  dbsi1e.o dbesi1.o  \
	dbsk0e.o  dbesk0.o  dbsk1e.o  dbesk1.o  \
	d9b0mp.o  d9b1mp.o  \
	dbesj0.o  dbesj1.o  dbesjn.o \
	dbesy0.o  dbesy1.o  dbesyn.o \
	d9aimp.o  daie.o dai.o  dbie.o  dbi.o  \
	ei.o  \
	jacobip.o  chebyshevt.o  laguerrel.o  hermiteh.o  \
	legendrep.o  legendrepa.o  harmonicy.o  \
	powi.o  powii.o  \
	clnrel.o  \
	cgamma.o  clngam.o  c9lgmc.o  \
	threejm.o  threejj.o  clebsch.o  sixj.o 

# implicit rule uses CXXFLAGS for compilation
CXXFLAGS=$(CCOPTS)

ALL: ${FCNOBJ} OBJLIST

${FCNOBJ}: ${FCNINC}

test: test.cc ${FCNOBJ}
	${CCC} -I ${MATPACK}/../include $(CCOPTS) -o $@ $@.cc ${FCNOBJ} ${LINKLIBS}

###########################################################################
### make a file "Objects" with a list of all ".o" files
###########################################################################

OBJLIST:
	@echo ${FCNOBJ} > Objects

###########################################################################
### exceptions
###########################################################################

###########################################################################
### make clean 								###
###########################################################################
 
clean:
	rm -f *.o core a.out

###########################################################################
### make protected save subdirectory
###########################################################################

save:
	@mkdir -p SAVE
	@chmod 755 SAVE
	@find *.c *.cc *.h *.doc *.tbl *.scn *.xpm *.cmp *.tex *.sty *.html \
	 *.xmgr *.xvgr Makefile* makefile* *.ps README COPYRIGHT > Files ; \
	 for i in `cat Files` ; do \
	     gzip -c $$i > SAVE/$$i.gz ; \
	 done 
	@chmod 500 SAVE

###########################################################################
