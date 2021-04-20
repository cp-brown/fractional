# Locations of source code, object files, and module files
SRC = ./src
OBJS = ./objs
MODS = ./mods

# Fortran compiler
F90 = ifort
FLAGS = -i8 -warn errors -fpp -fpic -heap-arrays -nogen-interfaces -module $(MODS)
#DEBUG = -Wall -Wextra -Wimplicit-interface -fPIC -fmax-errors=1 -g -fcheck=all -fbacktrace

OBJS = $(OBJS)/utils.o $(OBJS)/integration.o $(OBJS)/integration_trap_real.o $(OBJS)/integration_trap_cmplx.o
OBJS0D = $(OBJS)/inversion_0d.o $(OBJS)/invert1_0d.o $(OBJS)/invert2_0d.o $(OBJS)/invert3_0d.o
OBJSFEM1D = $(OBJS)/sparse_utils.o $(OBJS)/fem.o $(OBJS)/fem_mat_1d.o $(OBJS)/fem_vec_1d.o
OBJSFEM2D = ${OBJSFEM1D} $(OBJS)/fem_vec_2d.o $(OBJS)/fem_mat_2d.o
OBJS1D = $(OBJS)/inversion_1d.o $(OBJS)/invert1_1d.o $(OBJS)/invert2_1d.o $(OBJS)/invert3_1d.o
OBJS2D = $(OBJS)/inversion_2d.o $(OBJS)/invert1_2d.o $(OBJS)/invert2_2d.o $(OBJS)/invert3_2d.o

# For linking with MKL and BLAS95
MKLROOT = /opt/intel/oneapi/mkl/latest
MKLINCL = -I$(MKLROOT)/include/intel64/ilp64
MKLLINK = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl
BLAS95 = $(MKLROOT)/lib/intel64/libmkl_blas95_ilp64.a



TARGETS = 2009ex1 2003ex3 2009ex2
.PHONY: all clean
all: $(TARGETS)

clean:
	rm -f $(OBJS)/*.o $(MODS)/*.mod $(MODS)/*.smod

2009ex1: 2009ex1.f90 ${OBJS} ${OBJS0D}
	$(F90) $(FLAGS) $(MKLINCL) -o $@ $^ $(BLAS95) $(MKLLINK) -shared
	mv 2009ex1 ./2009ex1/

2003ex3: 2003ex3.f90 $(OBJS) $(OBJSFEM1D) $(OBJS1D)
	$(F90) $(FLAGS) $(MKLINCL) -o $@ $^ $(BLAS95) $(MKLLINK) -shared
	mv 2003ex3 ./2003ex3/

2009ex2: 2009ex2.f90 $(OBJS) $(OBJSFEM2D) $(OBJS2D)
	$(F90) $(FLAGS) $(MKLINCL) -o $@ $^ $(BLAS95) $(MKLLINK) -shared
	mv 2009ex2 ./2009ex2/

# ?????
%.o: %.f90
	$(F90) $(FLAGS) $(MKLINCL) -c -o $@ $<
