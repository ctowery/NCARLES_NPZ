# MODIFY THESE AS YOU NEED/LIKE ----------------------------------------------------------
SRC_DIR := ../code
MPI := mpif90
# note the use of = and not := here, allows deferring definitions of COMMON etc. to later
FFLAGS = $(COMMON) $(OPT2) # from most hardcore debugging to hardcore optimization, choose DBG3,2,1 or OPT1,2,3,4

# These flags are specific to gfortran ----------------------------------------------------
# I like verbose flags in a makefile, because why not when you aren't typing them out yourself?
# For gfortran:
# -r8 => -fdefault-real-8 -fdefault-double-8, this sets both `real` and `double precision` to 64 bits
# -132 => -ffixed-line-length-132, this lets all fixed-format *.f files have 132-character lines (default is 72)
COMMON := -fdefault-real-8 -fdefault-double-8 -ffixed-line-length-132 -ffree-line-length-none -std=legacy -fPIC -pipe
LDFLAGS := -lm

# Og allows some O1 optimizations while making debugging cleaner/clearer than O0
DBG1 := -Og -g1 -ffpe-trap=invalid,zero,overflow -fmax-errors=10 
DBG2 := -Og -g2 -ffpe-trap=invalid,zero,overflow,underflow -Wall -fcheck=all,no-array-temps # -g2 == -g
DBG3 := -O0 -g3 -ffpe-trap=invalid,zero,overflow,underflow -Wall -fcheck=all

OPT1 := -O1
OPT2 := -O2
OPT3 := -O3
OPT4 := -O3 -ffast-math -fno-protect-parens

# ----------------------------------------------------------------------------------------
SRCS := $(shell find $(SRC_DIR) -name '*.f' -or -name '*.f90')
INC_DIRS := $(shell find $(SRC_DIR) -type d)
vpath %.f90 $(INC_DIRS)  # necessary to compile the modules without explicitly specifying their directories

MOD_SRCS := inputs.f90 pars.f90 con_data.f90 con_stats.f90 fftwk.f90 fields.f90 tracerbc.f90 reaction.f90
MODS := $(MOD_SRCS:.f90=.mod)

# The final make step. Call with just `make` or explicitly as `make lespmi`.
lesmpi : $(MODS) 
	$(MPI) $(FFLAGS) $(SRCS) -o $@ $(LDFLAGS)

# make step for module files
%.mod : %.f90
	$(MPI) $(FFLAGS) -c $<


.PHONY : clean realclean cleansrc

# this just removes things made by the `make lesmpi` command in your current directory
clean: 
	rm -rf *.dSYM lesmpi

# this should remove everything in your current directory
realclean:
	rm -rf *.o *.mod *.dSYM lesmpi

# this removes things in the source directory. If compiling correctly these shouldn't exist,
# but if messing around with the makefile, you might accidentally compile something "in-place"
# within the source directory
cleansrc:
	rm -rf $(shell find $(SRC_DIR) -name '*.o' -or -name '*.mod' -or -name '*.dSYM')