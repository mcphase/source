#*********************************************************************
#
# File:          Makefile
#
# Project:       McPhase
#
# Description:   Makefile for the program package
#
# CHECKLIST - how to make a new distribution of McPhase:
#   -> edit always in mcphase~mcphase, 
#                                 make git
#                                 make
#               
#               IF NECESSARY: *********************************
#               .  oputdated: get newest version from repository (in mcphase~mcphase):
#                  git pull origin  master6_1
#               .  compile mcphaseexplorer 
#                              mcphaseexplorer compiled using netbeans if modified, in order 
#                               to use it: click on project right click package as zip. look in 
#                               Output/Netbeans_mpe/dist/ and unzip mcphaseexplorer.zip
#                               - store current bin/mcphaseexplorer as backup
#                               - unzip it 
#                               - move ions directory from old version to new and replace old version
#                               in mcphas/bin/mcphaseexplorer* with new version.
#                               - to avoid question on reload of .-jpg use own compiled  org-netbeans-modules-image.jar file:
#                               cd bin ; mv mcphaseexplorer mcphaseexplorer_old ; unzip ../Output/Netbeansproject_mpe/McPhaseExplorer/dist/mcphaseexplorer.zip
#                               mv mcphaseexplorer_old/ions  mcphaseexplorer/ ; cp ../Output/Netbeansproject_mpe/org-netbeans-modules-image.jar mcphaseexplorer/ide/modules/org-netbeans-modules-image.jar
#
#               .   make java [ compile java in Output/mcphase~mcphase/bin/jar, and make git to install it 
#                              (usually not required because java needs recompile only when changes are made]
#               .  update manual and latex in mcphase~mcphase/doc and create manual.pdf / html version with latex2html:
#                               make git ; make manual
#			        note: on MAC latex2html needs edit vi /opt/homebrew/Cellar/latex2html/2025/l2hconf.pm
#                                               and insert 
#$LATEX = '/Library/TeX/texbin/pdflatex';        # LaTeX $PDFLATEX = '/Library/TeX/texbin/pdflatex';     # pdfLaTeX
#$LUALATEX = '/Library/TeX/texbin/luatex';       # LuaLaTeX $DVILUALATEX = '/Library/TeX/texbin/dviluatex'; # dviLuaLaTeX
#$DVIPS = '/Library/TeX/texbin/dvips';   # dvips $DVIPNG = '/Library/TeX/texbin/dvipng'; # dvipng
#$PDFTOCAIRO = '/opt/homebrew/bin/pdftocairo';   # pdf to svg converter $PS2PDF = '/opt/homebrew/bin/ps2pdf';   # ps to pdf converter
#$PDFCROP = '/Library/TeX/texbin/pdfcrop';       # pdfcrop $GS = '/opt/homebrew/opt/ghostscript/bin/gs';   # GhostScript
#                                       on MAC latex2html:  if preview.sty is missing ---> "sudo tlmgr install preview"
#                                   --> this will create online manual in directory manual
#                                       copy this manual to the webpage: mkdir manual  |   put -r manual     (puts all files recursively)
#               .   outdated: push all changes to github, in mcphase~mcphase:
#                   git add bin/src/formfactor.c (to add new source files in case this is needed)
#                   git commit -a -m 'some changes in charges.c singleion.c'
#                   git push  -u origin master6_1
#               OBLIGATORY: *********************************************
#               .  update version information in version, version.tex
#               .  update innosetup (version information) // izpack.xl
#		.  type:  make  package
#                               or make windows (for windows) and make tgz (for linux)
#	           (does automatically make clean, recompile, creates in $HOME McPhaseSetup6_1.exe and McPhaseSetup6_1.tgz 
#                   and converts all *.pl files dos2unix in bin, perhaps in examples,
#                    -k keep going even if errors occur)
#               TESTING: *************************************************
#                     run demo in linux and windows
#                      - test demo 
#                      - in demo: mctest test*.bat
#               update WEB PAGES: ****************************************
#                       - check all information and add new, update list of papers on mcphas, references
#                       - put .tgz and .exe and if available additional software on web
#                       - test if mcphase can be installed on different computer systems
#               set GIT MASTER ******************************************* 
#                        ... after the next commit we should create a new branch by
#                        git checkout -b master6_2
#                         outdated: git push origin master6_2
#
#               Note on compiling for windows on McOS: make sure MacPorts is installed from https://www.macports.org/install.php 
#                                       and update it with : sudo port -v selfupdate
#                                       using this the windows gnu compilers can be installed on Mac: 
#                                                            sudo port install mingw-w64
#
# 	    gfortran library not found: try to find it with: find / -name "*libgfortran*" 2>/dev/null
#				and then do: 
#				export GFORTRANLIB=-L/opt/homebrew/Cellar/gcc/16.2.0/lib/gcc/current/
#           gfortran lSystem not found - then do:
#                               export LIBRARY_PATH="$LIBRARY_PATH:/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib"
#
#	        Note on compiling with dos:    copy mcphas (except Output)
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
#
#    Screen Saver
#
#It's pretty easy to make your own screen saver in IrfanView, really, all it is is a slideshow - sorry, 
# no transitional effects or background music for this. If you have recommendations for more 
# advanced software with these types of features, make a post in our forum. For clarity's sake, 
# I'm using IrfanView 3.85 - I would expect the process to be quite similar in recent past, 
# and near future versions of this popular program.
#
#Choose the image files with which you want to create a screen saver.
#Open IrfanView from the "Start" menu.
#In IrfanView, from the Menu Bar choose "File", then "SlideShow".
#Select image files and add them to your slideshow - they can be from anywhere on your computer. 
# Put them in the order that you want them to display. Unless you decide to use the random image 
# order options, then it would not matter in which order your files are added to the list.
# Once you've added all your files, click the "Save as EXE/SCR file" button. For a screen saver, 
# you want to choose "Create SCR file" in the next dialog box.
#Then select "Create" and the screen saver should be magically created for you.
# Author(s):     M. Rotter
#
#  Last Update:	  20.11.2026
#
#**********************************************************************
# On macOS, use Homebrew GCC's current library location unless supplied externally.
ifeq ($(shell uname -s),Darwin)
  ifeq ($(origin GFORTRANLIB),undefined)
    BREW_GCC_PREFIX := $(shell brew --prefix gcc 2>/dev/null)
    ifneq ($(strip $(BREW_GCC_PREFIX)),)
      GFORTRANLIB := -L$(BREW_GCC_PREFIX)/lib/gcc/current
    endif
  endif
endif
export GFORTRANLIB

include bin/src/Makefile.common

gitdir = Output/mcphase~mcphase

# Extract the numeric version from a line such as:
#   #define MCPHASVERSION "mcphas 6.1"
VERSION := $(shell sed -nE \
	's/.*MCPHASVERSION[[:space:]]+"[^"]*([0-9]+\.[0-9]+)".*/\1/p' \
	"version" | head -n 1)
VERSION_TAG := $(subst .,_,$(VERSION))

.PHONY: all allwin git manual vector vectorwin package windows tgz \
        clean cleanexe java winno help

.PHONY: help
help:
	@echo "McPhase targets:"
	@echo "  make              Build the main programs and examples"
	@echo "  make allwin       Cross-compile Windows programs"
	@echo "  make manual       Build the manual"
	@echo "  make package      Build Windows and Linux packages"
	@echo "  make clean        Remove generated build files"
	@echo "  make git          get newer files using rsync from repository in Output/mcphase~mpchase/"


all: vector larpack functions cfield mcphase ic1ion phonon examples tutorial bfk bcfph cowan 
allwin: vectorwin larpackwin functionswin cfieldwin mcphasewin ic1ionwin phononwin exampleswin tutorial bfkwin bcfphwin cowanwin


git: $(gitdir)/*
	cp Makefile Makefile.sav
#	cp -ru $(gitdir)/* ./
	rsync -avu --filter='exclude .git'  $(gitdir)/ ./
	cp Makefile.sav Makefile

manual:    
	cd doc && $(MAKE)

vector: 
	cd bin/src/vector && $(MAKE)

vectorwin: 
	cd bin/src/vector && $(MAKE) cross64=1

larpack:
	cd bin/ic1ion_module/libs ; $(MAKE)

larpackwin:
	cd bin/ic1ion_module/libs ; $(MAKE) cross64=1

functions: 
	cd bin/src/functions && $(MAKE)

functionswin: 
	cd bin/src/functions && $(MAKE) cross64=1

cfield: vector 
	cd bin/cf1ion_module && $(MAKE)

cfieldwin: vectorwin 
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

mcphase: vector cfield 
	cd bin/src && $(MAKE)

mcphasewin: vectorwin cfieldwin 
	cd bin/src && $(MAKE) cross64=1

examples  : vector 
	cd ./examples ; make

exampleswin  : vectorwin 
	cd ./examples ; $(MAKE) cross64=1

tutorial : vector
	cd ./tutorial ; make

java    : 
	cd ./bin/jar ; make

package :  
	$(MAKE) manual windows  tgz  

windows: 
	$(MAKE) manual clean cleanexe allwin 
	$(MAKE) clean 
	unix2dos ./bin/*.pl ./demo/*.bat
	substitute sleep wait ./demo/*.bat 
	unix2dos ./examples/cecu2a/fit/watch*.bat
	unix2dos ./examples/tmnic2/calc.bat
	unix2dos ./examples/coo/calc.bat
	unix2dos ./examples/dycu2iwata/calc.bat
	unix2dos ./examples/gd3gao6/calc.bat
	unix2dos ./examples/helix_spinwave/calc.bat
	unix2dos ./examples/ho2ti2o7/*.bat
	unix2dos ./examples/la2coo4/calc.bat
	unix2dos ./examples/Ce3p_chain_cfphonon/calc.bat
	unix2dos ./examples/Ce3p_tetragonalprim_cfphonon/calc.bat
	unix2dos ./examples/CeAl2_cfphonon_cfstrict/phonons/calc.bat
	unix2dos ./examples/CeAl2_cfphonon_cfstrict/CF/calc.bat
	unix2dos ./examples/Cu2OSeO3/4interacting_tetrahedrons_noDM/calc.bat
	unix2dos ./examples/Cu2OSeO3/nocluster_spinwave_model_noDM/calc.bat
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
	unix2dos ./examples/prni2si2/calc.bat
	unix2dos ./examples/pupd3/calc.bat
	unix2dos ./examples/Ru3p_create_sipf/calc.bat ./examples/Ru3p_create_sipf/calcsta.bat
	unix2dos ./examples/testic1ion/test.bat
	unix2dos ./examples/tmcu2_cf_phonon/calc.bat
	unix2dos ./examples/tungsten_phonons/calc.bat
	unix2dos ./examples/upd3/calc.bat
	unix2dos ./tutorial/07documentation_logbooks/calc.bat
#	/Applications/IzPack/bin/compile izpack.xl -o $(HOME)/windows.jar
#	python /Applications/IzPack/utils/wrappers/izpack2exe/izpack2exe.py --file=$(HOME)/windows.jar --no-upx --with-jdk=bin/zulu26.30.11-ca-jdk26.0.1-win_x64 --output=$(HOME)/mcph.exe
#	rm $(HOME)/windows.jar
	$(MAKE) winno

winno:
	dot_clean -mv ./
#	docker run --rm -i -v "$$PWD:/work" amake/innosetup innosetup_mac.iss
	cd bin/iss ; $(MAKE) clean ; $(MAKE) ; $(MAKE) clean
	mv Output/McPhaseSetup*.exe $(HOME)/

tgz : 
	$(MAKE) manual cleanexe clean all 
	$(MAKE) clean 
	dos2unix ./bin/*.pl
	dos2unix ./demo/*.bat ./demo/demo
	substitute wait sleep ./demo/*.bat 
	dos2unix ./examples/cecu2a/fit/watch*.bat
	dos2unix ./examples/tmnic2/calc.bat
	dos2unix ./examples/coo/calc.bat
	dos2unix ./examples/dycu2iwata/calc.bat
	dos2unix ./examples/gd3gao6/calc.bat
	dos2unix ./examples/helix_spinwave/calc.bat
	dos2unix ./examples/ho2ti2o7/*.bat
	dos2unix ./examples/la2coo4/calc.bat
	dos2unix ./examples/Ce3p_chain_cfphonon/calc.bat
	dos2unix ./examples/Ce3p_tetragonalprim_cfphonon/calc.bat
	dos2unix ./examples/CeAl2_cfphonon_cfstrict/phonons/calc.bat
	dos2unix ./examples/CeAl2_cfphonon_cfstrict/CF/calc.bat
	dos2unix ./examples/Cu2OSeO3/4interacting_tetrahedrons_noDM/calc.bat
	dos2unix ./examples/Cu2OSeO3/nocluster_spinwave_model_noDM/calc.bat
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
	dos2unix ./examples/prni2si2/calc.bat
	dos2unix ./examples/pupd3/calc.bat
	dos2unix ./examples/Ru3p_create_sipf/calc.bat ./examples/Ru3p_create_sipf/calcsta.bat
	dos2unix ./examples/testic1ion/test*.bat
	dos2unix ./examples/tmcu2_cf_phonon/calc.bat
	dos2unix ./examples/tungsten_phonons/calc.bat
	dos2unix ./examples/upd3/calc.bat
	dos2unix ./tutorial/07documentation_logbooks/calc.bat
	dot_clean -mv ./
	cd ../;tar --exclude='./mcphas/bin/zulu*' \
		--exclude=mcphas/bin/Perl* \
		--exclude=mcphas/Output* --exclude=mcphas/bin/*.exe \
		-cvf $(HOME)/mcph.tar mcphas/* \
		;tar -rvf $(HOME)/mcph.tar mcphas/bin/zulu27.30.109-ca-jdk27.0.1.0.3-macosx_aarch64-runtime/* \
		;tar -rvf $(HOME)/mcph.tar mcphas/bin/zulu25.38.109-ca-jdk25.0.5.0.3-linux_x64-runtime/* \
		;cd ./mcphas  
	cd $(HOME);gzip mcph.tar;mv mcph.tar.gz McPhaseSetup$(VERSION_TAG).tgz

clean:
	rm -f ./Makefile.sav
	cd ./doc ; $(MAKE) clean
	cd ./tutorial ; $(MAKE) clean
	cd ./examples ; $(MAKE) clean
	cd ./bin/jar ; $(MAKE) clean
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
	rm -vf bin/cond bin/clusterize bin/charges bin/coq2jjj \
		bin/jjj2j bin/mcdispit bin/mcphasit bin/anisotropyit bin/singleion bin/cfield \
		bin/spins  bin/spinsfromq  bin/sinlgeion \
                bin/densplt bin/pointc bin/formfactor\
                bin/radwavfunc bin/addj bin/reduce_unitcell \
                bin/mcdiff   bin/mf2fe bin/fediff bin/inimenu \
                bin/cf1ion_module/so1ion.so \
                bin/ic1ion  bin/icf1ion bin/so1ion \
                bin/ic1ion_module/ic1ion.so \
                bin/bfk \
                bin/RCN2K bin/RCN36K bin/RCG11K \
                bin/bcfph 
                