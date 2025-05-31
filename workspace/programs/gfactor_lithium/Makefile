FC = gfortran
FFLAGS = -Wall -O2
FLIBS = -lbspgf -lbspqed -lslatec -llapack

HOMDIR = /home/sriman/Documents/Thesis
WRKDIR = $(HOMDIR)/workspace
LIBDIR = $(WRKDIR)/include/lib
MODDIR = $(WRKDIR)/include/mod
PRJDIR = $(WRKDIR)/programs/gfactor_lithium
BINDIR = $(PRJDIR)/bin
SRCDIR = $(PRJDIR)/src

g_delbrueck: \
	$(SRCDIR)/g_delbrueck.f90
	$(FC) $(FFLAGS) -L$(LIBDIR) -I$(MODDIR) -o $(BINDIR)/$@ $< $(FLIBS)

g_dirac: \
	$(SRCDIR)/g_dirac.f90
	$(FC) $(FFLAGS) -L$(LIBDIR) -I$(MODDIR) -o $(BINDIR)/$@ $< $(FLIBS)

g_propagator_loop: \
	$(SRCDIR)/g_propagator_loop.f90
	$(FC) $(FFLAGS) -L$(LIBDIR) -I$(MODDIR) -o $(BINDIR)/$@ $< $(FLIBS)

g_uehling: \
	$(SRCDIR)/g_uehling.f90
	$(FC) $(FFLAGS) -L$(LIBDIR) -I$(MODDIR) -o $(BINDIR)/$@ $< $(FLIBS)

g_uehling_all: \
	$(SRCDIR)/g_uehling_all.f90
	$(FC) $(FFLAGS) -L$(LIBDIR) -I$(MODDIR) -o $(BINDIR)/$@ $< $(FLIBS)

g_wichmann_kroll: \
	$(SRCDIR)/g_wichmann_kroll.f90
	$(FC) $(FFLAGS) -L$(LIBDIR) -I$(MODDIR) -o $(BINDIR)/$@ $< $(FLIBS)
