#Unix makefile for fortran-file	

# put the name of the target program here
TARGET = brush.weak.sphere.dna # the list of source files
SRC =    precision.f90 mathconst.f90 physconst.f90 globals.f90 molecule.f90 chains.f90 volume.f90 rands.f90 myutils.f90 parameter.f90  field.f90 fenergy.f90 initcha.f90  myio.f90 cadenas.f90 fcnDNA.f90  rota.sphere.f90  init.f90 chaingenerator.f90 kinsolsolver.f90  solver.f90 L2norm.f90 brush.weak.sphere.csrange.f90
# fenergy.f90 after field.f90
# some definitions
SHELL = /bin/bash
# get git version
GIT_VERSION := $(shell git describe --abbrev=6 --dirty --always --tags)

is_quest = 

ifeq ($(shell hostname),gadol)

FFLAGS= -cpp -DVERSION=\"$(GIT_VERSION)\" -fbounds-check -Warray-bounds -Wall #-O3

LDFLAGS=-L/opt/local/sundials-2.5.0-openmpi/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/lib/gcc/x86_64-linux-gnu/4.6 -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../.. -lgfortran -lm -lgcc_s -lquadmath

LFFLAGS=$(LDFLAGS)

FF= gfortran

else ifeq ($(shell hostname),pear)

FFLAGS= -O3 -cpp -DVERSION=\"$(GIT_VERSION)\"
LDFLAGS=-L/opt/local/sundials-2.5.0-openmpi/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/lib/gcc/x86_64-linux-gnu/4.6 -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../.. -lgfortran -lm -lgcc_s -lquadmath

LFFLAGS=$(LDFLAGS)
FF= gfortran


else ifeq ($(shell hostname),master)

FFLAGS= -O3 -cpp -DVERSION=\"$(GIT_VERSION)\"
LDFLAGS=-L/shared/software/sundials-2.5.0-openmpi/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/lib/gcc/x86_64-linux-gnu/4.6 -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../.. -lgfortran -lm -lgcc_s -lquadmath

LFFLAGS=$(LDFLAGS)
FF= gfortran

else ifeq ($(shell hostname),chiquita)


FFLAGS= -O3 -cpp -DVERSION=\"$(GIT_VERSION)\"
LDFLAGS=-L/opt/local/sundials-2.5.0-openmpi-atlas/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/lib/gcc/x86_64-linux-gnu/4.6.1 -L/usr/lib/gcc/x86_64-linux-gnu/4.6.1/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.6.1/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.6.1/../../.. -lgfortran -lm -lgcc_s -lquadmath

LFFLAGS=$(LDFLAGS)
FF= gfortran

else ifeq ($(shell hostname),orange)

FFLAGS= -O3 -cpp -DVERSION=\"$(GIT_VERSION)\"

LDFLAGS=-L/opt/local/sundials-2.5.0-openmpi-atlas/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/lib/gcc/x86_64-linux-gnu/4.6.1 -L/usr/lib/gcc/x86_64-linux-gnu/4.6.1/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.6.1/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.6.1/../../.. -lgfortran -lm -lgcc_s -lquadmath

LFFLAGS=$(LDFLAGS)
FF= gfortran





else ifeq ($(shell hostname),master.bw01.bme.northwestern.edu)

FFLAGS= -O3 -cpp -DVERSION=\"$(GIT_VERSION)\"

LDFLAGS= -L/export/apps/sundials-2.5.0-openmpi/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/lib/gcc/x86_64-redhat-linux/4.4.6 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../.. -lgfortranbegin -lgfortran -lm -lgcc_s

LFFLAGS=$(LDFLAGS)

FF= gfortran

else ifeq ($(shell hostname),quser13)

	is_quest = yes

else ifeq ($(shell hostname),quser10)

	is_quest = yes

else 


FFLAGS= -fbounds-check -Warray-bounds -cpp -DVERSION=\"$(GIT_VERSION)\"

LDFLAGS=  -L/opt/local/kinsol-2.6.0/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/local/Cellar/gcc/4.9.2_1/lib/gcc/4.9/gcc/x86_64-apple-darwin13.4.0/4.9.2 -L/usr/local/Cellar/gcc/4.9.2_1/lib/gcc/4.9/gcc/x86_64-apple-darwin13.4.0/4.9.2/../../.. -lgfortran -lquadmath -lm


LFFLAGS=$(LDFLAGS)

FF= gfortran

endif

ifdef is_quest 


FFLAGS= -O3 -cpp -DVERSION=\"$(GIT_VERSION)\"

LDFLAGS=  -L/home/rna878/kinsol/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/opt/intel/composerxe-2011.3.174/tbb/lib/intel64/cc4.1.0_libc2.4_kernel2.6.16.21 -L/opt/intel/composerxe-2011.3.174/mkl/lib/intel64 -L/opt/intel/composerxe-2011.3.174/ipp/lib/intel64 -L/opt/intel/composerxe-2011.3.174/compiler/lib/intel64 -L/hpc/opt/intel/composerxe-2011.3.174/compiler/lib/intel64 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7/../../.. -L/lib64 -L/lib -L/usr/lib64 -L/usr/lib -limf -lm -lifport -lifcore -lsvml -lipgo -lirc -lpthread -lirc_s -ldl


LFFLAGS=$(LDFLAGS)

FF= ifort

endif



all:	$(TARGET)

$(TARGET): $(SRC:.f90=.o)
	$(FF) -o $(TARGET) $(SRC:.f90=.o)  $(LDFLAGS) $(LFLAGS)

#$(SRC:.f90=.o): 	
#	${FF} -c ${FFLAGS}  $(SRC) 

%.o: %.f90
#.f90.o :
	${FF} ${FFLAGS} -c $(SRC)

install: all
	cp $(TARGET) ~/bin


clean:	
	@rm -f $(SRC:.f90=.o) $(SRC:.f90=.d) $(TARGET) *~ *.mod

realclean: clean
	@rm -f .depend

depend dep:
	@$(FF)  $(CFLAGS) -MM $(SRC) > .depend 

ifeq (.depend, $(wildcard .depend))
include .depend
endif


















































