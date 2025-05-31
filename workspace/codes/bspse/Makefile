FC = gfortran
FFLAGS = -Wall -O2

HOMDIR = /home/sriman/Documents/Thesis
WRKDIR = $(HOMDIR)/workspace
LIBDIR = $(WRKDIR)/include/lib
MODDIR = $(WRKDIR)/include/mod
PRJDIR = $(WRKDIR)/codes/bspse
OBJDIR = $(PRJDIR)/obj
SRCDIR = $(PRJDIR)/src

all: \
	$(OBJDIR)/bspse_manypotential.o \
	$(OBJDIR)/bspse_onepotential.o \
	$(OBJDIR)/bspse_zeropotential.o

lib:
	ar rcs $(LIBDIR)/libbspse.a $(OBJDIR)/*.o

bspse_manypotential: $(OBJDIR)/bspse_manypotential.o
$(OBJDIR)/bspse_manypotential.o: \
	$(SRCDIR)/bspse_manypotential.f90
	$(FC) $(FFLAGS) -J$(MODDIR) -o $@ -c $<

bspse_onepotential: $(OBJDIR)/bspse_onepotential.o
$(OBJDIR)/bspse_onepotential.o: \
	$(SRCDIR)/bspse_onepotential.f90
	$(FC) $(FFLAGS) -J$(MODDIR) -o $@ -c $<

bspse_zeropotential: $(OBJDIR)/bspse_zeropotential.o
$(OBJDIR)/bspse_zeropotential.o: \
	$(SRCDIR)/bspse_zeropotential.f90
	$(FC) $(FFLAGS) -J$(MODDIR) -o $@ -c $<
