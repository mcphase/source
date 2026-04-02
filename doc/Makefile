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

all: manual

version.tex: manual.tex ../version 
	sed -e 's/#define.*VERSION *"//;s/"//' ../version > version.tex


manual.aux: manual.tex version.tex $(wildcard *.bib) $(wildcard *.tex) 
	latex manual.tex
	makeindex manual.idx

manual.bbl: $(wildcard ./*.bib) manual.tex  $(wildcard ./*.tex)
	latex manual
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


manual.pdf: manual.tex version.tex manual.bbl $(wildcard ./*.bib) $(wildcard ./*.tex)
	pdflatex  manual 
	pdflatex  manual 

manual: manual.pdf
	rm -rf ./manual/*
	latex2html -local_icons manual.tex 

clean :
	rm -vf manual.ps manual.dvi manual.aux manual.log manual.bbl manual.blg manual.toc version.tex manual.idx manual.ind 



