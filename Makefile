# Locations of source code, object files, and module files
SRC = ./src
OBJS = ./objs
MODS = ./mods

# Fortran compiler
F90 = ifort
DEBUG = -g -traceback -gen-interfaces -warn interfaces -check
FLAGS = -i8 -warn errors -fpp -heap-arrays -module $(MODS) $(DEBUG)

OBJSGEN = $(OBJS)/utils.o $(OBJS)/integration.o $(OBJS)/integration_trap_real.o $(OBJS)/integration_trap_cmplx.o
OBJS0D = $(OBJS)/inversion_0d.o $(OBJS)/invert1_0d.o $(OBJS)/invert2_0d.o $(OBJS)/invert3_0d.o
OBJSFEM1D = $(OBJS)/sparse_utils.o $(OBJS)/fem.o $(OBJS)/fem_mat_1d.o $(OBJS)/fem_vec_1d.o
OBJSFEM2D = ${OBJSFEM1D} $(OBJS)/fem_vec_2d.o $(OBJS)/fem_mat_2d.o
OBJS1D = $(OBJS)/inversion_1d.o $(OBJS)/invert1_1d.o $(OBJS)/invert2_1d.o $(OBJS)/invert3_1d.o
OBJS2D = $(OBJS)/inversion_2d.o $(OBJS)/invert1_2d.o $(OBJS)/invert2_2d.o $(OBJS)/invert3_2d.o
OBJSMKL = $(OBJS)/mkl_spblas.o $(OBJS)/mkl_dss.o

# For linking with MKL and BLAS95
MKLROOT = /opt/intel/oneapi/mkl/latest
MKLINCL = -I$(MKLROOT)/include/intel64/ilp64
MKLLINK = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl
BLAS95 = $(MKLROOT)/lib/intel64/libmkl_blas95_ilp64.a



TARGETS = 2009ex1_exec 2003ex3_exec 2009ex2_exec
.PHONY: all clean
all: $(TARGETS)

clean:
	rm -f $(OBJS)/*.o $(MODS)/*.mod $(MODS)/*.smod

2009ex1_exec: ${SRC}/2009ex1.f90 ${OBJSGEN} ${OBJS0D}
	$(F90) $(FLAGS) -o $@ $^
	mv 2009ex1_exec ./2009ex1/

2003ex3_exec: ${SRC}/2003ex3.f90 $(OBJSGEN) $(OBJSMKL) $(OBJSFEM1D) $(OBJS1D)
	$(F90) $(FLAGS) $(MKLINCL) -o $@ $^ $(BLAS95) $(MKLLINK)
	mv 2003ex3_exec ./2003ex3/

2009ex2_exec: ${SRC}/2009ex2.f90 $(OBJSGEN) $(OBJSMKL) $(OBJSFEM2D) $(OBJS2D)
	$(F90) $(FLAGS) $(MKLINCL) -o $@ $^ $(BLAS95) $(MKLLINK)
	mv 2009ex2_exec ./2009ex2/

$(OBJS)/%.o: $(SRC)/%.f90
	$(F90) $(FLAGS) $(MKLINCL) -c -o $@ $<

$(OBJS)/mkl_spblas.o: $(MKLROOT)/include/mkl_spblas.f90
	$(F90) $(FLAGS) $(MKLINCL) -c -o $@ $<

$(OBJS)/mkl_dss.o: $(MKLROOT)/include/mkl_dss.f90
	$(F90) $(FLAGS) $(MKLINCL) -c -o $@ $<