#*********************************************************************
#
# File:          Makefile
#
# Project:       McPhase
#
# Description:   Makefile for the program package
#
# Checklist - how to make a new distribution of McPhase:
#               .  compile mcphaseexplorer under windows and test
#               .  get newest version from repository
#                  git pull ssh://martin_rotter@git.kenai.com/mcphase~mcphase master
#               .  update version information in version, version.tex
#               .  update innosetup (version information)
#               .  update manual and latex -> manual.pdf
#               .   -> edit always in mcphase~mcphase, then
#                                 make git
#                                 make
#               .   push all changes to Kenai, e.g.
#                   git add bin/src/formfactor.c (to add new source files in case this is needed)
#                   git commit -a -m 'some changes in charges.c singleion.c'
#                   git push ssh://martin_rotter@git.kenai.com/mcphase~mcphase master5_2
#               . under linux: type make clean
#                              make unreleased_remove
#               .              make java [usually not required because java needs recompile only when changes are made]
#                              cd doc and delete directory manual
#                              latex2html -local_icons manual.tex  
#                                         [maybe necessary to adapt manual.tex style, usepackage
#                                          so that latex can run without problems on it]
#                                   --> this will create online manual in directory manual
#                                       copy this manual to the webpage: mkdir manual  |   put -r manual     (puts all files recursively)
#		.              compile linux (see install: make all)
#		.              make clean
#               . linux testing:
#                      - test demo and use stored examples to compare outputfiles using 
#                        perl important_scripts\compare_mcphase_results.pl examplesold examplesnew
#                      - (test examples) - test testic1ion
#		. linux create mcphas*_*.tgz:
#		.        make -k package
#	           (does automatically make unreleased_remove clean and tgz and converts all *.pl files dos2unix 
#                   in bin, perhaps in examples and unreleased_restore, -k keep going even if errors occur)
#               .  rename linux distribution file  Output/mcph.tar.gz to e.g. mcph5_2.tgz
#	        . change to dos:   make unreleased_remove and copy mcphas (except Output)
#                                   to c:\msys64/home/rotter/  
#                                    mingw win32 shell  (from startup menu) 
#                                    cd mcphas
#                                    . lin.bat                                                                     
#                                    make clean 
#                                    make 
#                                    make clean
#                                   create win-setup mcphX_Y.exe using
#                                   inno setup compiler: innosetup.iss in this directory
#                                     (rename file setup found in mcphas/output/setup.exe into mcphX_Y.exe)
#                                  - install windows&linux and test demo + examples
#                                 make unreleased_restore
#               . update WEB PAGES: - check all information and add new, update list of papers on mcphas, references
#                                   - put .tgz and .exe and if available additional software on web
#                                   - test if mcphase can be installed on different computer systems
#               . set GIT MASETER  ... after the next commit we should create a new branch by
#                              git checkout -b master5_3
#                              git push ssh://martin_rotter@git.kenai.com/mcphase~mcphase master5_3
#
#    Screen Saver
#
#It's pretty easy to make your own screen saver in IrfanView, really, all it is is a slideshow - sorry, no transitional effects or background music for this. If you have recommendations for more advanced software with these types of features, make a post in our forum. For clarity's sake, I'm using IrfanView 3.85 - I would expect the process to be quite similar in recent past, and near future versions of this popular program.
#
#Choose the image files with which you want to create a screen saver.
#Open IrfanView from the "Start" menu.
#In IrfanView, from the Menu Bar choose "File", then "SlideShow".
#Select image files and add them to your slideshow - they can be from anywhere on your computer. Put them in the order that you want them to display. Unless you decide to use the random image order options, then it would not matter in which order your files are added to the list.
#Once you've added all your files, click the "Save as EXE/SCR file" button. For a screen saver, you want to choose "Create SCR file" in the next dialog box.
#Then select "Create" and the screen saver should be magically created for you.
# Author(s):     M. Rotter
#
#  Last Update:	  15.11.2015
#
#**********************************************************************

include bin/src/Makefile.common

gitdir = Output/mcphase~mcphase




all: vector functions cfield mcphase examples tutorial bfk bfkp

git: $(gitdir)/*
	cp Makefile Makefile.sav
#	cp -ru $(gitdir)/* ./
	rsync -avu $(gitdir)/ ./
	cp Makefile.sav Makefile

vector: 
	cd bin/src/vector && $(MAKE)

functions: 
	cd bin/src/functions && $(MAKE)

cfield: vector ic1ion
	cd bin/cf1ion_module && $(MAKE)

ic1ion: vector 
	cd bin/ic1ion_module && $(MAKE)

phonon: vector 
	cd bin/phonon_module && $(MAKE)

bfk:   
	cd bin/bfk_src && $(MAKE)

bfkp:   
	cd bin/bfkp_src && $(MAKE)

mcphase: vector cfield phonon 
	cd bin/src && $(MAKE)

examples  : vector 
	cd ./examples ; make

tutorial : vector
	cd ./tutorial ; make

java    : 
	cd ./bin/jar ; make

unreleased_remove : 
	mv -n ./bin/src/clusterize.c ./Output/clusterize_notreleased.c
	cp -n ./Output/clusterize_pleasefund.c ./bin/src/clusterize.c
	rm -f ./bin/clusterize*

unreleased_restore : 
	mv ./Output/clusterize_notreleased.c ./bin/src/clusterize.c 

package : unreleased_remove all clean tgz unreleased_restore

tgz :
	dos2unix ./bin/*.pl
	dos2unix ./demo/*.bat ./demo/demo
	dos2unix ./examples/cecu2a/fit/watch*.bat
	dos2unix ./examples/coo/calc.bat
	dos2unix ./examples/dycu2iwata/calc.bat
	dos2unix ./examples/gd3gao6/calc.bat
	dos2unix ./examples/helix_spinwave/calc.bat
	dos2unix ./examples/ho2ti2o7/calc.bat
	dos2unix ./examples/la2coo4/calc.bat
	dos2unix ./examples/LaCoO3_podlesnyak_polaron/calc.bat
	dos2unix ./examples/LaCoO3_podlesnyak_polaron/2ions/calc.bat
	dos2unix ./examples/LaCoO3_podlesnyak_polaron/3ions/calc.bat
	dos2unix ./examples/LaCoO3_podlesnyak_polaron/5ions/calc.bat
	dos2unix ./examples/LaCoO3_podlesnyak_polaron/7ions/calc.bat
	dos2unix ./examples/lumno3/calc.bat
	dos2unix ./examples/ndba2cu3o7/calcsta.bat
	dos2unix ./examples/NiO/calc.bat
	dos2unix ./examples/Pr3Pd20Si6/calc.bat
	dos2unix ./examples/prni2b2c/fit/watch*.bat ./examples/prni2b2c/fit/calcsta 
	dos2unix ./examples/prni2b2c/powder_magnon.bat
	dos2unix ./examples/prni2si2/calc.bat
	dos2unix ./examples/pupd3/calc.bat
	dos2unix ./examples/testic1ion/test*.bat
	dos2unix ./examples/tungsten_phonons/calc.bat
	dos2unix ./examples/upd3/calc.bat
	dos2unix ./tutorial/07documentation_logbooks/calc.bat
	cd ../;tar --exclude=mcphas/bin/jre1.8.0_121/* --exclude=mcphas/bin/Perl* \
		--exclude=mcphas/Output* --exclude=mcphas/bin/*.exe \
		-cvf ./mcphas/Output/mcph.tar mcphas/* \
		;cd ./mcphas  
	cd ./Output;gzip mcph.tar
	unix2dos ./bin/*.pl ./demo/*.bat
	unix2dos ./examples/cecu2a/fit/watch*.bat
	unix2dos ./examples/coo/calc.bat
	unix2dos ./examples/dycu2iwata/calc.bat
	unix2dos ./examples/gd3gao6/calc.bat
	unix2dos ./examples/helix_spinwave/calc.bat
	unix2dos ./examples/ho2ti2o7/calc.bat
	unix2dos ./examples/la2coo4/calc.bat
	unix2dos ./examples/LaCoO3_podlesnyak_polaron/calc.bat
	unix2dos ./examples/LaCoO3_podlesnyak_polaron/2ions/calc.bat
	unix2dos ./examples/LaCoO3_podlesnyak_polaron/3ions/calc.bat
	unix2dos ./examples/LaCoO3_podlesnyak_polaron/5ions/calc.bat
	unix2dos ./examples/LaCoO3_podlesnyak_polaron/7ions/calc.bat
	unix2dos ./examples/lumno3/calc.bat
	unix2dos ./examples/ndba2cu3o7/calcsta.bat
	unix2dos ./examples/NiO/calc.bat
	unix2dos ./examples/Pr3Pd20Si6/calc.bat
	unix2dos ./examples/prni2b2c/fit/watch*.bat ./examples/prni2b2c/fit/calcsta 
	unix2dos ./examples/prni2b2c/powder_magnon.bat
	unix2dos ./examples/prni2si2/calc.bat
	unix2dos ./examples/pupd3/calc.bat
	unix2dos ./examples/testic1ion/test.bat
	unix2dos ./examples/tungsten_phonons/calc.bat
	unix2dos ./examples/upd3/calc.bat
	unix2dos ./tutorial/07documentation_logbooks/calc.bat

clean:
	rm -f ./Makefile.sav
	cd ./tutorial ; make clean
	cd ./examples ; make clean
	cd ./bin/jar ; make clean
	cd ./bin/src/vector && $(MAKE) clean
	cd ./bin/src/functions && $(MAKE) clean
	cd ./bin/phonon_module && $(MAKE) clean
	cd ./bin/cf1ion_module && $(MAKE) clean
	cd ./bin/ic1ion_module && $(MAKE) cleanall
	cd ./bin/src && $(MAKE) clean

cleanexe:
	rm -vf bin/addj.exe bin/charges.exe bin/coq2jjj.exe \
		bin/mcdispit.exe bin/singleion.exe bin/cfield.exe \
		bin/cond.exe bin/jjj2j.exe bin/mcphasit.exe bin/spins.exe \
                bin/chrgplt.exe bin/pointc.exe bin/spinsfromq.exe \
                bin/densplt.exe \
                bin/mcdiff.exe  \
                bin/ic1ion.exe bin/so1ion.exe \
                bin/fediff.exe bin/mf2fe.exe \
                bin/formfactor.exe bin/radwavfunc.exe bin/clusterize.exe
	rm -vf bin/addj bin/charges bin/coq2jjj \
		bin/mcdispit bin/singleion bin/cfield \
		bin/cond bin/jjj2j bin/mcphasit bin/spins \
                bin/chrgplt bin/pointc bin/spinsfromq \
                bin/densplt \
                bin/mcdiff bin/cf1ion_module/cfield.so \
                bin/ic1ion bin/so1ion \
                bin/ic1ion_module/ic1ion.so \
                bin/fediff bin/mf2fe \
                bin/formfactor bin/radwavfunc bin/clusterize
                