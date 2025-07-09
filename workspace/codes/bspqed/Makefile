FC = gfortran
FFLAGS = -Wall -O2

HOMDIR = /home/hrishikeshsriman-narayan/Documents/Masters/Master-Thesis--Weak_Int
WRKDIR = $(HOMDIR)/workspace
LIBDIR = $(WRKDIR)/include/lib
MODDIR = $(WRKDIR)/include/mod
PRJDIR = $(WRKDIR)/codes/bspqed
OBJDIR = $(PRJDIR)/obj
SRCDIR = $(PRJDIR)/src

all: \
	$(OBJDIR)/bspqed_coefficients.o \
	$(OBJDIR)/bspqed_delbrueck.o \
	$(OBJDIR)/bspqed_discretization.o \
	$(OBJDIR)/bspqed_dkb.o \
	$(OBJDIR)/bspqed_fourier_density.o \
	$(OBJDIR)/bspqed_function.o \
	$(OBJDIR)/bspqed_integrals.o \
	$(OBJDIR)/bspqed_momentum.o \
	$(OBJDIR)/bspqed_physconsts.o \
	$(OBJDIR)/bspqed_potentials.o \
	$(OBJDIR)/bspqed_quadrature.o \
	$(OBJDIR)/bspqed_uehling.o \
	$(OBJDIR)/bspqed_utilities.o \
	$(OBJDIR)/bspqed_wichmann_kroll.o

lib:
	ar rcs $(LIBDIR)/libbspqed.a $(OBJDIR)/*.o

bspqed_coefficients: $(OBJDIR)/bspqed_coefficients.o
$(OBJDIR)/bspqed_coefficients.o: \
	$(SRCDIR)/bspqed_coefficients.f90
	$(FC) $(FFLAGS) -J$(MODDIR) -o $@ -c $<

bspqed_delbrueck: $(OBJDIR)/bspqed_delbrueck.o
$(OBJDIR)/bspqed_delbrueck.o : \
	$(SRCDIR)/bspqed_delbrueck.f90 \
	$(OBJDIR)/bspqed_function.o \
	$(OBJDIR)/bspqed_quadrature.o
	$(FC) $(FFLAGS) -J$(MODDIR) -o $@ -c $<

bspqed_discretization: $(OBJDIR)/bspqed_discretization.o
$(OBJDIR)/bspqed_discretization.o: \
	$(SRCDIR)/bspqed_discretization.f90 \
	$(OBJDIR)/bspqed_dkb.o \
	$(OBJDIR)/bspqed_function.o
	$(FC) $(FFLAGS) -J$(MODDIR) -o $@ -c $<

bspqed_dkb: $(OBJDIR)/bspqed_dkb.o
$(OBJDIR)/bspqed_dkb.o: \
	$(SRCDIR)/bspqed_dkb.f90 \
	$(OBJDIR)/bspqed_function.o
	$(FC) $(FFLAGS) -J$(MODDIR) -o $@ -c $<

bspqed_fourier_density: $(OBJDIR)/bspqed_fourier_density.o
$(OBJDIR)/bspqed_fourier_density.o: \
	$(SRCDIR)/bspqed_fourier_density.f90 \
	$(OBJDIR)/bspqed_function.o
	$(FC) $(FFLAGS) -J$(MODDIR) -o $@ -c $<

bspqed_function: $(OBJDIR)/bspqed_function.o
$(OBJDIR)/bspqed_function.o: \
	$(SRCDIR)/bspqed_function.f90
	$(FC) $(FFLAGS) -J$(MODDIR) -o $@ -c $<

bspqed_integrals: $(OBJDIR)/bspqed_integrals.o
$(OBJDIR)/bspqed_integrals.o : \
	$(SRCDIR)/bspqed_integrals.f90 \
	$(OBJDIR)/bspqed_dkb.o
	$(FC) $(FFLAGS) -J$(MODDIR) -o $@ -c $<

bspqed_momentum: $(OBJDIR)/bspqed_momentum.o
$(OBJDIR)/bspqed_momentum.o: \
	$(SRCDIR)/bspqed_momentum.f90
	$(FC) $(FFLAGS) -J$(MODDIR) -o $@ -c $<

bspqed_physconsts: $(OBJDIR)/bspqed_physconsts.o
$(OBJDIR)/bspqed_physconsts.o: \
	$(SRCDIR)/bspqed_physconsts.f90
	$(FC) $(FFLAGS) -J$(MODDIR) -o $@ -c $<

bspqed_potentials: $(OBJDIR)/bspqed_potentials.o
$(OBJDIR)/bspqed_potentials.o: \
	$(SRCDIR)/bspqed_potentials.f90 \
	$(OBJDIR)/bspqed_function.o \
	$(OBJDIR)/bspqed_uehling.o
	$(FC) $(FFLAGS) -J$(MODDIR) -o $@ -c $<

bspqed_quadrature: $(OBJDIR)/bspqed_quadrature.o
$(OBJDIR)/bspqed_quadrature.o: \
	$(SRCDIR)/bspqed_quadrature.f90
	$(FC) $(FFLAGS) -J$(MODDIR) -o $@ -c $<

bspqed_uehling: $(OBJDIR)/bspqed_uehling.o
$(OBJDIR)/bspqed_uehling.o : \
	$(SRCDIR)/bspqed_uehling.f90 \
	$(OBJDIR)/bspqed_function.o
	$(FC) $(FFLAGS) -J$(MODDIR) -o $@ -c $<

bspqed_utilities: $(OBJDIR)/bspqed_utilities.o
$(OBJDIR)/bspqed_utilities.o : \
	$(SRCDIR)/bspqed_utilities.f90 \
	$(OBJDIR)/bspqed_function.o \
	$(OBJDIR)/bspqed_integrals.o
	$(FC) $(FFLAGS) -J$(MODDIR) -o $@ -c $<

bspqed_wichmann_kroll: $(OBJDIR)/bspqed_wichmann_kroll.o
$(OBJDIR)/bspqed_wichmann_kroll.o : \
	$(SRCDIR)/bspqed_wichmann_kroll.f90 \
	$(OBJDIR)/bspqed_function.o
	$(FC) $(FFLAGS) -J$(MODDIR) -o $@ -c $<
