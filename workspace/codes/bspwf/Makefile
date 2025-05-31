FC = gfortran
FFLAGS = -Wall -O2

HOMDIR = /home/sriman/Documents/Thesis
WRKDIR = $(HOMDIR)/workspace
LIBDIR = $(WRKDIR)/include/lib
MODDIR = $(WRKDIR)/include/mod
PRJDIR = $(WRKDIR)/codes/bspwf
OBJDIR = $(PRJDIR)/obj
SRCDIR = $(PRJDIR)/src

all: \
	$(OBJDIR)/bspwf_constants.o \
	$(OBJDIR)/bspwf_dirac.o \
	$(OBJDIR)/bspwf_function.o \
	$(OBJDIR)/bspwf_hartree_fock.o \
	$(OBJDIR)/bspwf_knots.o \
	$(OBJDIR)/bspwf_potentials.o \
	$(OBJDIR)/bspwf_quadrature.o \
	$(OBJDIR)/bspwf_radial_integrals.o \
	$(OBJDIR)/bspwf_schroedinger.o \
	$(OBJDIR)/bspwf_utilities.o

lib:
	ar rcs $(LIBDIR)/libbspwf.a $(OBJDIR)/*.o

bspwf_constants: $(OBJDIR)/bspwf_constants.o
$(OBJDIR)/bspwf_constants.o: \
	$(SRCDIR)/bspwf_constants.f90
	$(FC) $(FFLAGS) -J$(MODDIR) -o $@ -c $<

bspwf_dirac: $(OBJDIR)/bspwf_dirac.o
$(OBJDIR)/bspwf_dirac.o: \
	$(SRCDIR)/bspwf_dirac.f90 \
	$(OBJDIR)/bspwf_constants.o \
	$(OBJDIR)/bspwf_function.o
	$(FC) $(FFLAGS) -J$(MODDIR) -o $@ -c $<

bspwf_function: $(OBJDIR)/bspwf_function.o
$(OBJDIR)/bspwf_function.o: \
	$(SRCDIR)/bspwf_function.f90
	$(FC) $(FFLAGS) -J$(MODDIR) -o $@ -c $<

bspwf_hartree_fock: $(OBJDIR)/bspwf_hartree_fock.o
$(OBJDIR)/bspwf_hartree_fock.o: \
	$(SRCDIR)/bspwf_hartree_fock.f90 \
	$(OBJDIR)/bspwf_constants.o
	$(FC) $(FFLAGS) -J$(MODDIR) -o $@ -c $<

bspwf_knots: $(OBJDIR)/bspwf_knots.o
$(OBJDIR)/bspwf_knots.o: \
	$(SRCDIR)/bspwf_knots.f90
	$(FC) $(FFLAGS) -J$(MODDIR) -o $@ -c $<

bspwf_potentials: $(OBJDIR)/bspwf_potentials.o
$(OBJDIR)/bspwf_potentials.o: \
	$(SRCDIR)/bspwf_potentials.f90 \
	$(OBJDIR)/bspwf_function.o
	$(FC) $(FFLAGS) -J$(MODDIR) -o $@ -c $<

bspwf_quadrature: $(OBJDIR)/bspwf_quadrature.o
$(OBJDIR)/bspwf_quadrature.o: \
	$(SRCDIR)/bspwf_quadrature.f90
	$(FC) $(FFLAGS) -J$(MODDIR) -o $@ -c $<

bspwf_radial_integrals: $(OBJDIR)/bspwf_radial_integrals.o
$(OBJDIR)/bspwf_radial_integrals.o: \
	$(SRCDIR)/bspwf_radial_integrals.f90 \
	$(OBJDIR)/bspwf_constants.o
	$(FC) $(FFLAGS) -J$(MODDIR) -o $@ -c $<

bspwf_schroedinger: $(OBJDIR)/bspwf_schroedinger.o
$(OBJDIR)/bspwf_schroedinger.o: \
	$(SRCDIR)/bspwf_schroedinger.f90 \
	$(OBJDIR)/bspwf_function.o
	$(FC) $(FFLAGS) -J$(MODDIR) -o $@ -c $<

bspwf_utilities: $(OBJDIR)/bspwf_utilities.o
$(OBJDIR)/bspwf_utilities.o: \
	$(SRCDIR)/bspwf_utilities.f90
	$(FC) $(FFLAGS) -J$(MODDIR) -o $@ -c $<
