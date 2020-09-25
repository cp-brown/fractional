FC90 = mpifort

DEBUG = -Wall -Wextra -Wimplicit-interface -fPIC -fmax-errors=1 -g -fcheck=all -fbacktrace

TARGETS = 2009ex1 2003ex3 2009ex2

OBJS = utils.o \
       integration.o integration_trap_real.o integration_trap_cmplx.o
OBJSMKL = blas.o mkl_spblas.o mkl_dss.o
OBJS0D = inversion_0d.o invert1_0d.o invert2_0d.o invert3_0d.o
OBJSFEM1D = sparse_utils.o fem.o fem_mat_1d.o fem_vec_1d.o
OBJSFEM2D = ${OBJSFEM1D} fem_vec_2d.o fem_mat_2d.o
OBJS1D = inversion_1d.o invert1_1d.o invert2_1d.o invert3_1d.o
OBJS2D = inversion_2d.o invert1_2d.o invert2_2d.o invert3_2d.o

MKLROOT = /home/cpbrown/intel/compilers_and_libraries_2020.1.217/linux/mkl
BLAS95ROOT = ${MKLROOT}/interfaces/blas95
MKL = ${BLAS95ROOT}/lib/intel64/libmkl_blas95_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl -m64 -I${MKLROOT}/include

.PHONY: all clean

all: $(TARGETS)

clean:
	rm *.o *.mod *.smod

2009ex1: 2009ex1.f90 ${OBJS} ${OBJS0D}
	$(FC90) -o 2009ex1 $^
	mv 2009ex1 ../2009ex1/

2003ex3: 2003ex3.f90 $(OBJS) $(OBJSMKL) $(OBJSFEM1D) $(OBJS1D)
	$(FC90) -o $@ $^ $(MKL)
	mv 2003ex3 ../2003ex3/

2009ex2: 2009ex2.f90 $(OBJS) $(OBJSMKL) $(OBJSFEM2D) $(OBJS2D)
	$(FC90) $(DEBUG) -o $@ $^ $(MKL)
	mv 2009ex2 ../2009ex2/

%.o: %.f90
	$(FC90) -c -o $@ $<
