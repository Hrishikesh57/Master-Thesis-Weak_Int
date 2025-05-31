FC = gfortran
FFLAGS = -Wall -O2
FLIBS = -lbspwf -lslatec -llapack

HOMDIR = /home/sriman/Documents/Thesis
WRKDIR = $(HOMDIR)/workspace
LIBDIR = $(WRKDIR)/include/lib
MODDIR = $(WRKDIR)/include/mod
PRJDIR = $(WRKDIR)/programs/hartree
BINDIR = $(PRJDIR)/bin
SRCDIR = $(PRJDIR)/src

dhfcb: \
	$(SRCDIR)/dhfcb.f90
	$(FC) $(FFLAGS) -L$(LIBDIR) -I$(MODDIR) -o $(BINDIR)/$@ $< $(FLIBS)

dhfc: \
	$(SRCDIR)/dhfc.f90
	$(FC) $(FFLAGS) -L$(LIBDIR) -I$(MODDIR) -o $(BINDIR)/$@ $< $(FLIBS)

shf: \
	$(SRCDIR)/shf.f90
	$(FC) $(FFLAGS) -L$(LIBDIR) -I$(MODDIR) -o $(BINDIR)/$@ $< $(FLIBS)
