#################################################################################
#  										#
# Makefile for the MatPack 1.0 C++ vector class libary      by B.M.Gammel, 1995	#
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
### include common rules
###########################################################################

MATPACK=..
include $(MATPACK)/Makefile.common

INCLUDE= ../../include


.SUFFIXES: .cc .o

.c.o: 
	${CC} -c ${COPTS} -I ${INCLUDE} -I../../ic1ion_module $<

.cc.o: 
	${CCC} -c ${CCOPTS} -I ${INCLUDE} -I../../ic1ion_module $<



###########################################################################
### make all objects
###########################################################################

VECINC=	${INCLUDE}/vector.h

VECOBJ= matpack.o sortdbl.o sort2dbl.o sort3dbl.o\
	ivector.o imatrix.o \
	dvector.o dmatrix.o dmatmul.o \
	cvector.o cmatrix.o cmatmul.o \
	dludcomp.o dsvdcomp.o dortho.o \
	cludcomp.o \
	tred.o imtql.o imtql2.o rseigen.o cheigen.o chgeigen.o  \
	matlapl.o moment.o matexph.o matsqth.o \
	magicsquare.o

ALL: ${VECOBJ} OBJLIST

${VECOBJ}: ${VECINC}

###########################################################################
### make a file "Objects" with a list of all ".o" files
###########################################################################

OBJLIST:
	@echo ${VECOBJ} > Objects


###########################################################################
### make clean 								###
###########################################################################
 
clean:
	rm -vf ${VECOBJ} 

###########################################################################
