FC = gfortran
FFLAGS = -Wall -O2
FLIBS = -lbspse -lbspqed -lslatec -llapack

OMPFLAGS = -fopenmp

HOMDIR = /home/sriman/Documents/Thesis
WRKDIR = $(HOMDIR)/workspace
LIBDIR = $(WRKDIR)/include/lib
MODDIR = $(WRKDIR)/include/mod
PRJDIR = $(WRKDIR)/programs/selfen
BINDIR = $(PRJDIR)/bin
SRCDIR = $(PRJDIR)/src

many: \
	$(SRCDIR)/many.f90
	$(FC) $(FFLAGS) $(OMPFLAGS) -L$(LIBDIR) -I$(MODDIR) -o $(BINDIR)/$@ $< $(FLIBS)

one: \
	$(SRCDIR)/one.f90
	$(FC) $(FFLAGS) -L$(LIBDIR) -I$(MODDIR) -o $(BINDIR)/$@ $< $(FLIBS)

zero: \
	$(SRCDIR)/zero.f90
	$(FC) $(FFLAGS) -L$(LIBDIR) -I$(MODDIR) -o $(BINDIR)/$@ $< $(FLIBS)
