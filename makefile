# Makefile

include make.inc

EXE_DIR := bin/
SRC_DIR := src/
OBJ_DIR := obj/

# program file name
PROG := $(EXE_DIR)ptr.out
#PROG := $(EXE_DIR)interp.out

# compiler
FXX := ${COMPILER_NAME}

# compiler option

FFLAGS:=${OPTION}
# FFLAGS += -O0 -CB -traceback -g -fpe0 -check uninit -check pointers # -warn unused
#FFLAGS += -diag-disable=10121
#FFLAGS += -module $(OBJ_DIR)
FFLAGS += -J$(OBJ_DIR)

# library link
LIBS := -lhdf5_hl_fortran -lhdf5_hl

# suffix rule
.SUFFIXES: .f90 .F90 .o .mod
%.o: %.mod

# source file
SRC:=\
input_parser.f90\
unit_mod.f90\
const_mod.f90\
module_eos.f90\
io_mod.f90\
simdata3D_mod.f90\
particle_data_mod.f90\
particle_set_mod.f90\
divide.f90\
module_rho_ye.f90\
module_restart_hdf.f90\
module_flux_sample.f90\
condition_ejecta.f90\
subroutine_analysis.f90\
main.f90\
evolution_particle_3D.f90\
set_ejecta_inside_3D_divide.f90\
set_ejecta_uniform.f90\
recursive_division.f90\
module_store.f90\
ascii.f90\
analysis.f90\
print_data.f90\
partial_output_hdf.f90

# SRC:=\
# input_parser.f90\
# unit_mod.f90\
# module_eos.f90\
# simdata3D_mod.f90\
# interp_data.f90\
# interpolation.f90

SRC := $(addprefix $(SRC_DIR), $(SRC))

OBJ_FILES := $(addprefix $(OBJ_DIR),$(notdir $(SRC:.f90=.o)))
MOD_FILES := $(addprefix $(OBJ_DIR),$(notdir $(MOD:.f90=.mod)))

.PHONY: all dirs

all: dirs $(PROG)

objs : dirs $(OBJ_FILES)

dirs : $(EXE_DIR) $(OBJ_DIR)


$(EXE_DIR):
	mkdir -p $(EXE_DIR)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(PROG): $(OBJ_FILES)
	$(FXX) $(FFLAGS) -o $@ $(OBJ_FILES) $(OMPFLAGS) $(LIBS)

$(OBJ_DIR)%.o: $(SRC_DIR)%.f90
	$(FXX) $(FFLAGS) $(OMPFLAGS) -c $< -o $@

$(OBJ_DIR)%.o: $(SRC_DIR)%.F90
	$(FXX) $(FFLAGS) $(OMPFLAGS) -c $< -o $@

$(OBJ_DIR)%.mod: $(SRC_DIR)%.f90
	$(FXX) $(FFLAGS) $(OMPFLAGS) -c $< -o $@

$(OBJ_DIR)%.mod: $(SRC_DIR)%.F90
	$(FXX) $(FFLAGS) $(OMPFLAGS) -c $< -o $@

# clean rule
.PHONY: clean
clean:
	$(RM) $(OBJ_FILES) $(MOD_FILES) fort* 
