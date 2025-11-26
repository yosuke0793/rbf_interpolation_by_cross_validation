#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Makefile for DIFFERENTIAL EVOLUTION
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Name of executable files
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
 EXENAM = rbf_interpolation
 TARGET = $(HOME)/bin/$(EXENAM)
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Compiler selection ( Required, Default:FC=ifort )
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# INTEL FORTRAN 90 ( ifort )
  FC = ifort
#
# GNU FORTRAN 90 ( gfortran )
# FC = gfortran
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Compiler options ( Optional )
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Debugging [0] -> Optimization [4]
# OPTFLAG = 0
# OPTFLAG = 1
# OPTFLAG = 2
# OPTFLAG = 3
OPTFLAG = 4
#
# Fast option
# OPTFLAG = 5
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Advanced options for optimization ( Optional )
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#MCNFLAG = nut
#MCNFLAG = tan
MCNFLAG = iruka
# MCNFLAG = same
#MCNFLAG = potato
#MNCFLAG = pumpkin
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Intel Math Kernel Library ( Required )
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
 MKLDIR = /opt/intel/oneapi/mkl/latest/lib/intel64/
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Linker Options ( Required )
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
ifeq ($(FC),ifort)
	# INTEL FORTRAN 90 COMPILER
	LDFLAGS  = -qmkl
else ifeq ($(FC),gfortran)
	# GNU FORTRAN 90 COMPILER
	LDFLAGS  = -Wl,--start-group
	LDFLAGS += $(MKLDIR)libmkl_intel_lp64.a
	LDFLAGS += $(MKLDIR)libmkl_intel_thread.a
	LDFLAGS += $(MKLDIR)libmkl_core.a
	LDFLAGS += -ldl
	LDFLAGS += -Wl,--end-group -liomp5 -lpthread
endif
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Compiler options ( GNU Make v3.81+ )
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
ifeq ($(FC),ifort)
	ifeq ($(OPTFLAG),0)
		FFLAGS  = -g -O0 -traceback
		FFLAGS += -check all -warn all
		FFLAGS += -fpe0 -ftrapuv -gen_interfaces
		FFLAGS += -fstack-security-check
		FFLAGS += -fstack-protector-all
		FFLAGS += -e90 -e03
	else ifeq ($(OPTFLAG),1)
		FFLAGS  = -O1 -traceback
	else ifeq ($(OPTFLAG),2)
		FFLAGS  = -O2 -static -traceback
	else ifeq ($(OPTFLAG),3)
		FFLAGS  = -O3 -traceback -qopenmp
	else ifeq ($(OPTFLAG),4)
		FFLAGS  = -O3 -ipo -static -qopenmp
	else ifeq ($(OPTFLAG),5)
		FFLAGS  = -O3 -ipo -no-prec-div -static
	endif
else ifeq ($(FC),gfortran)
	ifeq ($(OPTFLAG),0)
		FFLAGS  = -g -O0 -fbacktrace -fbounds-check -Wall
		FFLAGS += -ffpe-trap=invalid,zero,overflow,underflow,denormal
		FFLAGS += -pedantic -std=gnu
		FFLAGS += -O -Wuninitialized
	else ifeq ($(OPTFLAG),1)
		FFLAGS  = -O1
	else ifeq ($(OPTFLAG),2)
		FFLAGS  = -O2
	else ifeq ($(OPTFLAG),3)
		FFLAGS  = -O3
	else ifeq ($(OPTFLAG),4)
		FFLAGS  = -O3
	else ifeq ($(OPTFLAG),5)
		FFLAGS  = -Ofast
	endif
endif
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Compiler options ( Advanced )
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
ifeq ($(FC),ifort)
	ifeq ($(MCNFLAG),nut)
		FFLAGS += -axAVX -xAVX
	else ifeq ($(MCNFLAG),tan)
		FFLAGS += -axCORE-AVX2 -xCORE-AVX2
	else ifeq ($(MCNFLAG),iruka)
		FFLAGS += -axCORE-AVX512 -xCORE-AVX512
	else ifeq ($(MCNFLAG),same)
		FFLAGS += -march=core-avx2
	else ifeq ($(MCNFLAG),potato)
		FFLAGS += -axSSE4.2 -xSSE4.2
	else ifeq ($(MCNFLAG),pumpkin)
		FFLAGS += -axCORE-AVX2 -xCORE-AVX2
	endif
endif
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  User defined function
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
define GENMV
	@if ls *genmod* 1> /dev/null 2>&1; then mv *genmod* $(OBJDIR); fi
endef
#
define MODMV
	@if ls *.mod 1> /dev/null 2>&1; then mv *.mod $(OBJDIR); fi
endef
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  MPI Library selection ( Required )
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
MPICOMD = $(shell which mpirun)
MPIPATH = $(or $(findstring /opt/openmpi, $(MPICOMD)), $(findstring /usr/local/openmpi, $(MPICOMD)))
#
ifeq ($(MPIPATH),/opt/openmpi)
	# OPEN MPI
	MPIFLAG = 1
else ifeq ($(MPIPATH),/usr/local/openmpi)
	# OPEN MPI
	MPIFLAG = 2
else
	# INTEL MPI ( mpifort: ifort, mpif90: gfortran )
	MPIFLAG = 3
endif
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  MPI Compiler selection ( GNU Make v3.81+ )
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
ifeq ($(MPIFLAG),1)
	ifeq ($(FC),ifort)
		FC = mpifort -diag-disable=10448
	else ifeq ($(FC),gfortran)
		FC = env OMPI_FC=gfortran mpifort
	endif
else ifeq ($(MPIFLAG),2)
	ifeq ($(FC),ifort)
		FC = mpifort -diag-disable=10448
	else ifeq ($(FC),gfortran)
		FC = mpif90
	endif
else ifeq ($(MPIFLAG),3)
	ifeq ($(FC),ifort)
		FC = mpifort -diag-disable=10448
	else ifeq ($(FC),gfortran)
		FC = mpif90
	endif
endif
#
ifeq ($(FC),mpifort)
	ifeq ($(shell which $(FC) > /dev/null 2>&1; echo $$?), 1)
		FC = mpif90
	endif
endif
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Intermediate directory path
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
OBJDIR = ./.intermediate
FFLAGS += -I$(OBJDIR)
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  File storing folders
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
vpath %.f90 ./module
vpath %.f90 ./main00
vpath %.f90 ./utilts
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  File names
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
OBJ001  = commonval.o fileinfo.o 
OBJ002  = util_common.o util_math.o util_rbf.o util_de.o
OBJ003  = main00.o
#
OBJECTS  = $(OBJ001) $(OBJ002) $(OBJ003)
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Make Target
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
all: $(TARGET)
.PHONY: all
#
$(TARGET): $(patsubst %,$(OBJDIR)/%,$(OBJECTS))
	@echo " "
	@echo " Linking ... "
	@$(FC) $^ $(LDFLAGS) -o $@
	@$(GENMV)
	@$(MODMV)
	@echo " "
	@echo "============================================="
	@echo " '$(EXENAM)' make done.  "
	@echo "============================================="
	@echo " "
#
$(OBJDIR)/%.o : %.f90
ifeq ($(EXITFLAG),1)
	$(error MPI Gfortran environment is not set)
else
	@mkdir -p $(HOME)/bin
	@mkdir -p $(@D)
	$(FC) -c $(FFLAGS) $< -o $@
endif
#
clean:
	@$(RM) $(TARGET)
	@$(RM) *.o *.mod *_genmod*
	@$(RM) $(OBJDIR)/*.o $(OBJDIR)/*.mod $(OBJDIR)/*_genmod*
#
