# Makefile

include make.inc

EXE_DIR := bin/
OBJ_DIR := obj/

SRC_DIR := src/
SRC2_DIR := ../src/

PROG := $(EXE_DIR)a.out
ifndef PROG
  $(error MAKE_TARGET is not set well. Specify ptr or interp)
endif
# compiler
FXX := ${COMPILER_NAME}

# compiler option

FFLAGS:=${OPTION}
# FFLAGS += -O0 -CB -traceback -g -fpe0 -check uninit -check pointers # -warn unused
#FFLAGS += -diag-disable=10121
FFLAGS += -module $(OBJ_DIR)
#FFLAGS += -J$(OBJ_DIR)

# library link
LIBS := -lhdf5_hl_fortran -lhdf5_hl

# suffix rule
.SUFFIXES: .f90 .F90 .o .mod
%.o: %.mod

# source file
MOD2 := module_eos.f90
SRC2 :=

MOD1:= input_parser.f90 unit_mod.f90 module_weak_interaction.f90 module_etanu.f90 #module_EOS_sekig.f90


SRC1:= make_table.f90 test_weak_rate.f90 test_ye_equil.f90 sort_median.f90

MOD :=
MOD += $(addprefix $(SRC2_DIR), $(MOD2))
MOD += $(addprefix $(SRC_DIR), $(MOD1))

SRC :=
SRC += $(addprefix $(SRC_DIR), $(SRC1))
SRC += $(addprefix $(SRC2_DIR), $(SRC2))

SRC := $(MOD) $(SRC)

OBJ_FILES := $(addprefix $(OBJ_DIR),$(notdir $(SRC:.f90=.o)))
MOD_FILES := $(addprefix $(OBJ_DIR),$(notdir $(MOD:.f90=.mod)))

.PHONY: all dirs test

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


$(OBJ_DIR)%.o: $(SRC2_DIR)%.f90
	$(FXX) $(FFLAGS) $(OMPFLAGS) -c $< -o $@

$(OBJ_DIR)%.o: $(SRC2_DIR)%.F90
	$(FXX) $(FFLAGS) $(OMPFLAGS) -c $< -o $@

$(OBJ_DIR)%.mod: $(SRC2_DIR)%.f90
	$(FXX) $(FFLAGS) $(OMPFLAGS) -c $< -o $@

$(OBJ_DIR)%.mod: $(SRC2_DIR)%.F90
	$(FXX) $(FFLAGS) $(OMPFLAGS) -c $< -o $@

test:
	@echo $(OBJ_FILES)
	@echo $(SRC)
	@echo $(SRC2)

# clean rule
.PHONY: clean
clean:
	$(RM) $(OBJ_FILES) $(MOD_FILES) fort* 
