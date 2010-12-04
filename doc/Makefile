#*********************************************************************
#
# File:          Makefile
#
# Project:       McPhase Manual
#
# Description:   Makefile for the program
#
# Author(s):     M. Rotter
#
# Original:	 August 2001
#
#**********************************************************************

all: manual.pdf  

manual.aux: manual.tex $(wildcard *.bib) $(wildcard *.tex) ../bin/src/dos_fileman/readme.txt
	sed -e 's/#define.*VERSION *"//;s/"//' ../version > version.tex
	latex manual.tex
	makeindex manual.idx

manual.bbl: $(wildcard ./*.bib) manual.tex manual.aux
	bibtex manual


manual.dvi: manual.tex manual.aux manual.bbl manual.toc $(wildcard ./*.eps) $(wildcard ./figsrc/*.ps)
	latex manual.tex
	latex manual.tex
	
manual.ps: manual.dvi $(wildcard ./*.eps)  $(wildcard ./figsrc/*.ps)
	dvips manual.dvi

#dp.ps: manual.ps
#	pstops  "2:0L@.65(21cm,0)+1L@.65(21cm,14.85cm)" $< > dp.ps

#dev.ps: dp.ps
#	pstops  "2:-1" $<  > dev.ps

#dodd.ps: dp.ps
#	pstops  "2:0" $<  > dodd.ps


manual.pdf: manual.ps
	ps2pdf  $<

clean :
	rm -vf manual.aux manual.log manual.bbl manual.blg manual.toc version.tex manual.idx manual.ind 



