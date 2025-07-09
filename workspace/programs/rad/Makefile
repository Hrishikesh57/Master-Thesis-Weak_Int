FC = gfortran
FFLAGS = -Wall -O2
FLIBS = -lbspwf -lbspgf -lbspqed -lslatec -llapack

HOMDIR = /home/hrishikeshsriman-narayan/Documents/Masters/Master-Thesis--Weak_Int
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
mewf : \
	$(SRCDIR)/mewf.f90
	$(FC) $(FFLAGS) -L$(LIBDIR) -I$(MODDIR) -o $(BINDIR)/$@ $< $(FLIBS)

bspwf_h: \
	$(SRCDIR)/bspwf_h.f90
	$(FC) $(FFLAGS) -L$(LIBDIR) -I$(MODDIR) -o $(BINDIR)/$@ $< $(FLIBS)
