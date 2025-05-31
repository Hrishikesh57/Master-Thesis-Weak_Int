# FC = gfortran
# FFLAGS = -Wall -O2
# FLIBS = -lbspgf -lbspqed -lslatec -llapack

# HOMDIR = /home/sriman/Documents/Thesis
# WRKDIR = $(HOMDIR)/workspace
# LIBDIR = $(WRKDIR)/include/lib
# MODDIR = $(WRKDIR)/include/mod
# PRJDIR = $(WRKDIR)/codes/bspgf
# OBJDIR = $(PRJDIR)/obj
# SRCDIR = $(PRJDIR)/src
# BINDIR = $(PRJDIR)/bin

# all: \
# 	$(OBJDIR)/bspgf_hydrogen.o \
# 	$(OBJDIR)/bspgf_boron.o \
# 	$(OBJDIR)/bspgf_lithium.o \
# 	$(OBJDIR)/practice.o \
# 	$(OBJDIR)/theoretical_radial_function_1s_2s.o \
# 	$(OBJDIR)/radial_integrals.o 


# lib:
# 	ar rcs $(LIBDIR)/libbspgf.a $(OBJDIR)/*.o

# bspgf_hydrogen: $(OBJDIR)/bspgf_hydrogen.o
# $(OBJDIR)/bspgf_hydrogen.o: \
# 	$(SRCDIR)/bspgf_hydrogen.f90
# 	$(FC) $(FFLAGS) -J$(MODDIR) -o $@ -c $<


# bspgf_boron: $(OBJDIR)/bspgf_boron.o
# $(OBJDIR)/bspgf_boron.o: \
# 	$(SRCDIR)/bspgf_boron.f90
# 	$(FC) $(FFLAGS) -J$(MODDIR) -o $@ -c $<

# bspgf_lithium: $(OBJDIR)/bspgf_lithium.o
# $(OBJDIR)/bspgf_lithium.o: \
# 	$(SRCDIR)/bspgf_lithium.f90
# 	$(FC) $(FFLAGS) -J$(MODDIR) -o $@ -c $<

# bspgf_practice: $(OBJDIR)/practice.o
# $(OBJDIR)/practice.o: \
# 	$(SRCDIR)/practice.f90
# 	$(FC) $(FFLAGS) -J$(MODDIR) -o $@ -c $<

# bspgf_theoretical_radial_function_1s_2s: $(OBJDIR)/theoretical_radial_function_1s_2s.o
# $(OBJDIR)/theoretical_radial_function_1s_2s.o: \
# 	$(SRCDIR)/theoretical_radial_function_1s_2s.f90
# 	$(FC) $(FFLAGS) -J$(MODDIR) -o $@ -c $<

# bspgf_radial_integrals: $(OBJDIR)/radial_integrals.o
# $(OBJDIR)/radial_integrals.o: \
# 	$(SRCDIR)/radial_integrals.f90
# 	$(FC) $(FFLAGS) -J$(MODDIR) -o $@ -c $<


# # Compiler and Flags# Compiler and Flags
FC = gfortran
FFLAGS = -Wall -O2

# Directories
HOMDIR = /home/sriman/Documents/Thesis
WRKDIR = $(HOMDIR)/workspace
LIBDIR = $(WRKDIR)/include/lib
MODDIR = $(WRKDIR)/include/mod
PRJDIR = $(WRKDIR)/codes/bspgf
OBJDIR = $(PRJDIR)/obj
SRCDIR = $(PRJDIR)/src

# Source Files
SRC_HYDROGEN = $(SRCDIR)/bspgf_hydrogen.f90
SRC_BORON    = $(SRCDIR)/bspgf_boron.f90
SRC_LITHIUM  = $(SRCDIR)/bspgf_lithium.f90
SRC_PRACTICE = $(SRCDIR)/practice.f90
SRC_RADIAL   = $(SRCDIR)/theoretical_radial_function_1s_2s.f90
SRC_RADIAL_INTEGRALS = $(SRCDIR)/radial_integrals.f90  # Added source file

# Object Files
OBJ_HYDROGEN = $(OBJDIR)/bspgf_hydrogen.o
OBJ_BORON    = $(OBJDIR)/bspgf_boron.o
OBJ_LITHIUM  = $(OBJDIR)/bspgf_lithium.o
OBJ_PRACTICE = $(OBJDIR)/practice.o
OBJ_RADIAL   = $(OBJDIR)/theoretical_radial_function_1s_2s.o
OBJ_RADIAL_INTEGRALS = $(OBJDIR)/radial_integrals.o    # Added object file

# Executables
EXEC_PRACTICE = $(PRJDIR)/practice1
EXEC_RADIAL   = $(PRJDIR)/theoretical_radial_function_1s_2s
EXEC_RADIAL_INTEGRALS = $(PRJDIR)/radial_integrals    # Added executable

# Libraries
LIB_BSPGF    = $(LIBDIR)/libbspgf.a
LIB_OTHER    = -lbspqed -lslatec -llapack -lbspwf

# Default Target
.PHONY: all
all: lib practice theoretical_radial_function_1s_2s radial_integrals    # Added radial_integrals

# Build Object Files
.PHONY: objects
objects: $(OBJ_HYDROGEN) $(OBJ_BORON) $(OBJ_LITHIUM) $(OBJ_PRACTICE) $(OBJ_RADIAL) $(OBJ_RADIAL_INTEGRALS)   # Added OBJ_RADIAL_INTEGRALS

# Compile Hydrogen Module
$(OBJ_HYDROGEN): $(SRC_HYDROGEN)
	$(FC) $(FFLAGS) -J$(MODDIR) -c $< -o $@

# Compile Boron Module
$(OBJ_BORON): $(SRC_BORON)
	$(FC) $(FFLAGS) -J$(MODDIR) -c $< -o $@

# Compile Lithium Module
$(OBJ_LITHIUM): $(SRC_LITHIUM)
	$(FC) $(FFLAGS) -J$(MODDIR) -c $< -o $@

# Compile Practice Program
$(OBJ_PRACTICE): $(SRC_PRACTICE)
	$(FC) $(FFLAGS) -J$(MODDIR) -c $< -o $@

# Compile Theoretical Radial Function Program
$(OBJ_RADIAL): $(SRC_RADIAL)
	$(FC) $(FFLAGS) -J$(MODDIR) -c $< -o $@

# Compile Radial Integrals Program
$(OBJ_RADIAL_INTEGRALS): $(SRC_RADIAL_INTEGRALS)
	$(FC) $(FFLAGS) -J$(MODDIR) -c $< -o $@    # Compiling with module directory

# Create Static Library
.PHONY: lib
lib: $(LIB_BSPGF)

$(LIB_BSPGF): $(OBJ_HYDROGEN) $(OBJ_BORON) $(OBJ_LITHIUM)
	ar rcs $@ $^

# Link Practice Executable
.PHONY: practice
practice: $(EXEC_PRACTICE)

$(EXEC_PRACTICE): $(OBJ_PRACTICE) $(LIB_BSPGF)
	$(FC) $(FFLAGS) -I$(MODDIR) -o $@ $< -L$(LIBDIR) $(LIB_OTHER) -lbspgf

# Link Theoretical Radial Function Executable
.PHONY: theoretical_radial_function_1s_2s
theoretical_radial_function_1s_2s: $(EXEC_RADIAL)

$(EXEC_RADIAL): $(OBJ_RADIAL)
	$(FC) $(FFLAGS) -I$(MODDIR) -o $@ $<

# Link Radial Integrals Executable
.PHONY: radial_integrals
radial_integrals: $(EXEC_RADIAL_INTEGRALS)

$(EXEC_RADIAL_INTEGRALS): $(OBJ_RADIAL_INTEGRALS) $(LIB_BSPGF)
	$(FC) $(FFLAGS) -I$(MODDIR) -o $@ $< -L$(LIBDIR) $(LIB_OTHER) -lbspgf

# Run Targets
.PHONY: run_practice run_theoretical_radial_function_1s_2s run_radial_integrals

run_practice: practice
	$(EXEC_PRACTICE)

run_theoretical_radial_function_1s_2s: theoretical_radial_function_1s_2s
	$(EXEC_RADIAL)

run_radial_integrals: radial_integrals
	$(EXEC_RADIAL_INTEGRALS)

# Clean Target
.PHONY: clean
clean:
	rm -f $(OBJDIR)/*.o $(EXEC_PRACTICE) $(EXEC_RADIAL) $(EXEC_RADIAL_INTEGRALS) $(LIB_BSPGF)
