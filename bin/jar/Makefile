#*********************************************************************
#
# File:          Makefile 
#
# Project:       McPhase
#
# Description:   Makefile for the java module
#
# Author(s):     M. Rotter
#
# Original:	 November 2011
#
#**********************************************************************

all : pack 

OBJS = IniConfigurator.java IniFile.java IniItem.java IniKey.java display.java \
       displaycontour.java displaydsigma.java displaymag.java \
       displaytext.java FileListener.java FileMonitor.java searchspacestatus.java \
       simannfitstatus.java Viewer.java ToolTipManager.java

compile : $(OBJS)
	javac ./*.java
#	gcj ./*.java

pack    : compile
	jar cvf mcphas.jar ./*.class
#	gjar cvf mcphas.jar ./*.class

clean  :
	rm -vf ./*.class
