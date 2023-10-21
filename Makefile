#*********************************************************************
#
# File:          Makefile
#
# Project:       McPhase
#
# Description:   Makefile for the program package
#
# CHECKLIST - how to make a new distribution of McPhase:
#   -> edit always in mcphase~mcphase, then
#                                 make git
#                                 make
#               
#               IF NECESSARY: *********************************
#               .  get newest version from repository
#                  git pull ssh://martin_rotter@git.kenai.com/mcphase~mcphase master
#               .  compile mcphaseexplorer under windows and test
#                              mcphaseexplorer compiled using netbeans if modified, in order 
#                               to use it: click on project right click package as zip. look in 
#                               Output/Netbeans_mpe/dist/ and unzip mcphaseexplorer.zip
#                               unzip it 
#                               move ions directory from old version to new and replace old version
#                               in mcphas/bin/mcphaseexplorer* with new version.
#               .   make java [ compile java in bin/jar,usually not required because java needs recompile only when changes are made]
#               .  update manual and latex -> manual.pdf
#                              cd doc and delete directory manual
#                              latex2html -local_icons manual.tex  
#                                         [maybe necessary to adapt manual.tex style, usepackage
#                                          so that latex can run without problems on it]
#                                   --> this will create online manual in directory manual
#                                       copy this manual to the webpage: mkdir manual  |   put -r manual     (puts all files recursively)
#               .   outdated: push all changes to Kenai, e.g.
#                   git add bin/src/formfactor.c (to add new source files in case this is needed)
#                   git commit -a -m 'some changes in charges.c singleion.c'
#                   git push ssh://martin_rotter@git.kenai.com/mcphase~mcphase master5_2
#               OBLIGATORY: *********************************************
#               .  update version information in version, version.tex
#               .  update innosetup (version information) // izpack.xl
#		.      start docker and type:  make  package
#                               or make windows (for windows) and make tgz (for linux)
#	           (does automatically make unreleased_remove clean, mcph.exe and mcph.tgz 
#                   and converts all *.pl files dos2unix 
#                   in bin, perhaps in examples and unreleased_restore,
#                    -k keep going even if errors occur)
#               .  rename linux distribution file  $HOME/mcph.tgz to e.g. mcph5_2.tgz
#               .  rename windows distribution file  $HOME/mcph.exe to e.g. mcph5_2.exe
#               TESTING: *************************************************
#                     run demo in linux and windows
#                      - test demo and use stored examples to compare outputfiles using 
#                        perl important_scripts\compare_mcphase_results.pl examplesold examplesnew
#                      - (test examples) - test testic1ion
#               update WEB PAGES: ****************************************
#                                   - check all information and add new, update list of papers on mcphas, references
#                                   - put .tgz and .exe and if available additional software on web
#                                   - test if mcphase can be installed on different computer systems
#               set GIT MASTER ******************************************* 
#                              ... after the next commit we should create a new branch by
#                              git checkout -b master5_5
#                              git push ssh://martin_rotter@git.kenai.com/mcphase~mcphase master5_5
#
#
#
#               Note on compiling for windows on McOS: make sure MacPorts is installed from https://www.macports.org/install.php 
#                                       and update it with : sudo port -v selfupdate
#                                       using this the windows gnu compilers can be installed on Mac: 
#                                                            sudo port install mingw-w64
#
#
#                                       [ alternative: install izpack and izpack utils - izpack2exe
#                                       search careful for 7zS.sxf ... in https://github.com/github/ghfw-build-extra/blob/master/7-Zip/7zSD.sfx
#                                          and substitute the file in utils/wrapper/izpack2exe/ with that
#                                           correct bugs in izpack2exe.py (format of config.txt errors)]
#	    gfortran library not found: try to find it with: find / -name "*libgfortran*" 2>/dev/null
#				and then do: 
#				export GFORTRANLIB=-L/opt/homebrew/Cellar/gcc/13.1.0/lib/gcc/current/
#	        Note on compiling with dos:   make unreleased_remove and copy mcphas (except Output)
#                                   to c:\msys64/home/rotter/  
#                                    mingw win64 shell  (from startup menu) 
#                                    cd mcphas
#                                    . lin.bat 
#				   find . -name  '.DS_Store' -delete
#				   find . -name  '._*' -delete                                                                   
#                                    make clean 
#                                    make 
#                                    make clean
#                                   create win-setup mcphX_Y.exe using
#                                   inno setup compiler: innosetup.iss in this directory
#                                     (rename file setup found in mcphas/output/setup.exe into mcphX_Y.exe)
#                                  - install windows&linux and test demo + examples, attention: remove
#                                       HOME/appdat/roaming/.mcphaseexplorer directory before installation
#                                 make unreleased_restore
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
#  Last Update:	  22.02.2023
#
#**********************************************************************

include bin/src/Makefile.common

gitdir = Output/mcphase~mcphase




all: vector functions cfield mcphase phonon examples tutorial bfk bcfph cowan
allwin: vectorwin functionswin cfieldwin mcphasewin phononwin exampleswin tutorial bfkwin bcfphwin cowanwin


git: $(gitdir)/*
	cp Makefile Makefile.sav
#	cp -ru $(gitdir)/* ./
	rsync -avu --filter='exclude .git'  $(gitdir)/ ./
	cp Makefile.sav Makefile

vector: 
	cd bin/src/vector && $(MAKE)

vectorwin: 
	cd bin/src/vector && $(MAKE) cross64=1

functions: 
	cd bin/src/functions && $(MAKE)

functionswin: 
	cd bin/src/functions && $(MAKE) cross64=1

cfield: vector ic1ion
	cd bin/cf1ion_module && $(MAKE)

cfieldwin: vectorwin ic1ionwin
	cd bin/cf1ion_module && $(MAKE) cross64=1

ic1ion: vector 
	cd bin/ic1ion_module && $(MAKE)

ic1ionwin: vectorwin
	cd bin/ic1ion_module && $(MAKE) cross64=1

phonon: vector 
	cd bin/phonon_module && $(MAKE)

phononwin: vectorwin
	cd bin/phonon_module && $(MAKE) cross64=1

bfk:   
	cd bin/bfk_src && $(MAKE)

bfkwin:   
	cd bin/bfk_src && $(MAKE) cross64=1

bcfph:   mcphase
	cd bin/bcfph_src && $(MAKE)

bcfphwin:  mcphasewin
	cd bin/bcfph_src && $(MAKE) cross64=1

cowan:   
	cd bin/cowan && $(MAKE)

cowanwin:   
	cd bin/cowan && $(MAKE) cross64=1

mcphase: vector cfield phonon 
	cd bin/src && $(MAKE)

mcphasewin: vectorwin cfieldwin phononwin 
	cd bin/src && $(MAKE) cross64=1

examples  : vector 
	cd ./examples ; make

exampleswin  : vectorwin 
	cd ./examples ; make cross64=1

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

package :  
	make windows  tgz  

windows: 
	make unreleased_remove clean cleanexe allwin 
	make clean 
#	/Applications/IzPack/bin/compile izpack.xl -o $(HOME)/windows.jar
#	python /Applications/IzPack/utils/wrappers/izpack2exe/izpack2exe.py --file=$(HOME)/windows.jar --no-upx --with-jdk=bin/jre1.8.0_121 --output=$(HOME)/mcph.exe
#	rm $(HOME)/windows.jar
	dot_clean -mv ./
	docker run --rm -i -v "$$PWD:/work" amake/innosetup innosetup_mac.iss
	mv Output/mysetup.exe $(HOME)/mcph.exe
	make unreleased_restore

tgz : 
	make unreleased_remove cleanexe clean all 
	make clean 
	dos2unix ./bin/*.pl
	dos2unix ./demo/*.bat ./demo/demo
	dos2unix ./examples/cecu2a/fit/watch*.bat
	dos2unix ./examples/cecu2ge2_cf_phonon_int/phonons/calc.bat ./examples/cecu2ge2_cf_phonon_int/DMD_method/calc.bat
	dos2unix ./examples/coo/calc.bat
	dos2unix ./examples/dycu2iwata/calc.bat
	dos2unix ./examples/gd3gao6/calc.bat
	dos2unix ./examples/helix_spinwave/calc.bat
	dos2unix ./examples/ho2ti2o7/calc.bat
	dos2unix ./examples/la2coo4/calc.bat
	dos2unix ./examples/Ce3p_chain_cfphonon/calc.bat
	dos2unix ./examples/Ce3p_tetragonalprim_cfphonon/calc.bat
	dos2unix ./examples/CeAl2_cfphonon_cfstrict/phonons/calc.bat
	dos2unix ./examples/CeAl2_cfphonon_cfstrict/DMD_method/calc.bat
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
	dos2unix ./examples/Ru3p_create_sipf/calc.bat ./examples/Ru3p_create_sipf/calcsta.bat
	dos2unix ./examples/tb2ti2o7_cfphonon_simple/DMD_method/calc.bat
	dos2unix ./examples/tb2ti2o7_cfphonon_simple/phonons/calc.bat
	dos2unix ./examples/test_cluster/calc.bat
	dos2unix ./examples/testic1ion/test*.bat
	dos2unix ./examples/tmcu2_cf_phonon/calc.bat
	dos2unix ./examples/tungsten_phonons/calc.bat
	dos2unix ./examples/upd3/calc.bat
	dos2unix ./tutorial/07documentation_logbooks/calc.bat
	dot_clean -mv ./
	cd ../;tar --exclude=mcphas/bin/jre1.8.0_121/* --exclude=mcphas/bin/Perl* \
		--exclude=mcphas/Output* --exclude=mcphas/bin/*.exe \
		-cvf $(HOME)/mcph.tar mcphas/* \
		;cd ./mcphas  
	cd $(HOME);gzip mcph.tar;mv mcph.tar.gz mcph.tgz
	unix2dos ./bin/*.pl ./demo/*.bat
	unix2dos ./examples/cecu2a/fit/watch*.bat
	unix2dos ./examples/cecu2ge2_cf_phonon_int/phonons/calc.bat ./examples/cecu2ge2_cf_phonon_int/DMD_method/calc.bat
	unix2dos ./examples/coo/calc.bat
	unix2dos ./examples/dycu2iwata/calc.bat
	unix2dos ./examples/gd3gao6/calc.bat
	unix2dos ./examples/helix_spinwave/calc.bat
	unix2dos ./examples/ho2ti2o7/calc.bat
	unix2dos ./examples/la2coo4/calc.bat
	unix2dos ./examples/Ce3p_chain_cfphonon/calc.bat
	unix2dos ./examples/Ce3p_tetragonalprim_cfphonon/calc.bat
	unix2dos ./examples/CeAl2_cfphonon_cfstrict/phonons/calc.bat
	unix2dos ./examples/CeAl2_cfphonon_cfstrict/DMD_method/calc.bat
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
	unix2dos ./examples/Ru3p_create_sipf/calc.bat ./examples/Ru3p_create_sipf/calcsta.bat
	unix2dos ./examples/tb2ti2o7_cfphonon_simple/DMD_method/calc.bat
	unix2dos ./examples/tb2ti2o7_cfphonon_simple/phonons/calc.bat
	unix2dos ./examples/test_cluster/calc.bat
	unix2dos ./examples/testic1ion/test.bat
	unix2dos ./examples/tmcu2_cf_phonon/calc.bat
	unix2dos ./examples/tungsten_phonons/calc.bat
	unix2dos ./examples/upd3/calc.bat
	unix2dos ./tutorial/07documentation_logbooks/calc.bat
	make unreleased_restore

clean:
	rm -f ./Makefile.sav
	cd ./doc ; make clean
	cd ./tutorial ; make clean
	cd ./examples ; make clean
	cd ./bin/jar ; make clean
	cd ./bin/src/vector && $(MAKE) clean
	cd ./bin/src/functions && $(MAKE) clean
	cd ./bin/phonon_module && $(MAKE) clean
	cd ./bin/cf1ion_module && $(MAKE) clean
	cd ./bin/ic1ion_module && $(MAKE) cleanall
	cd ./bin/src && $(MAKE) clean
	cd ./bin/bfk_src && $(MAKE) clean
	cd ./bin/bcfph_src && $(MAKE) clean
	cd ./bin/cowan && $(MAKE) clean

cleanexe:
	rm -vf bin/addj.exe bin/anisotropyit.exe bin/charges.exe bin/coq2jjj.exe \
		bin/mcdispit.exe bin/singleion.exe bin/cfield.exe \
		bin/cond.exe bin/jjj2j.exe bin/mcphasit.exe bin/spins.exe \
                bin/chrgplt.exe bin/pointc.exe bin/spinsfromq.exe \
                bin/densplt.exe \
                bin/mcdiff.exe bin/reduce_unitcell.exe  \
                bin/ic1ion.exe bin/icf1ion.exe bin/so1ion.exe \
                bin/fediff.exe bin/mf2fe.exe \
                bin/formfactor.exe bin/radwavfunc.exe bin/clusterize.exe bin/bfk.exe \
                bin/RCN2K.exe bin/RCN36K.exe bin/RCG11K.exe \
                bin/bcfph.exe
	rm -vf bin/addj bin/charges bin/coq2jjj \
		bin/mcdispit bin/singleion bin/cfield \
		bin/cond bin/jjj2j bin/mcphasit bin/spins \
                bin/chrgplt bin/pointc bin/spinsfromq \
                bin/densplt bin/icf1ion \
                bin/mcdiff bin/reduce_unitcell bin/cf1ion_module/cfield.so \
                bin/ic1ion bin/so1ion \
                bin/ic1ion_module/ic1ion.so \
                bin/fediff bin/mf2fe \
                bin/formfactor bin/radwavfunc bin/clusterize bin/bfk \
                bin/RCN2K bin/RCN36K bin/RCG11K \
                bin/bcfph 
                