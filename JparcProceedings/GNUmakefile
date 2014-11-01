# $Id: GNUmakefile,v 1.2 2006/11/13 08:57:40 kuzniak Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := JparcProceedings
G4TARGET := $(name)
G4EXLIB := true

CPPFLAGS += -I/usr/include/qt4/Qt -I/usr/include/qt4

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

#-----------------------------------------------------------
# environment setting for ROOT output
#-----------------------------------------------------------

# enable ROOT output option
CPPFLAGS += -D__OUTPUT_ROOTFILE__

# add compile flags and extra libs for ROOT output
CPPFLAGS += $(shell $(ROOTSYS)/bin/root-config --cflags)
	ROOTLIBS = $(shell $(ROOTSYS)/bin/root-config --glibs)
	EXTRALIBS += $(ROOTLIBS)

#-----------------------------------------------------------


.PHONY: all
all: lib bin

include $(G4INSTALL)/config/architecture.gmk

include $(G4INSTALL)/config/binmake.gmk
