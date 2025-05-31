FC = gfortran
FFLAGS = -Wall -O2
FLIBS = -lbspgf -lbspqed -lslatec -llapack -lbspwf

HOMDIR = /home/sriman/Documents/Thesis
WRKDIR = $(HOMDIR)/workspace
LIBDIR = $(WRKDIR)/include/lib
MODDIR = $(WRKDIR)/include/mod
PRJDIR = $(WRKDIR)/programs/rad
BINDIR = $(PRJDIR)/bin
SRCDIR = $(PRJDIR)/src

something: \
	$(SRCDIR)/something.f90
	$(FC) $(FFLAGS) -L$(LIBDIR) -I$(MODDIR) -o $(BINDIR)/$@ $< $(FLIBS)
practise: \
	$(SRCDIR)/practise.f90
	$(FC) $(FFLAGS) -L$(LIBDIR) -I$(MODDIR) -o $(BINDIR)/$@ $< $(FLIBS)
theoretical_radial_functions_1s_2s: \
	$(SRCDIR)/theoretical_radial_functions_1s_2s.f90
	$(FC) $(FFLAGS) -L$(LIBDIR) -I$(MODDIR) -o $(BINDIR)/$@ $< $(FLIBS)
 