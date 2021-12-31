
FC := gfortran

PROGRAM := BtoJ.x
SRC     := ./src/
OBJ     := ./obj/
MOD     := ./mod/
BIN     := ./bin/
DEST    := /usr/local/bin

# For use of large stack memory, add '-Wl,-stack_size,0x40000000,-stack_addr,0xf0000000'
# For debug, replace '-O3' with '-O0 -g -Wall -fbounds-check'
# When using LAPACK in MacOS, add following flag to FFLAGS (NOTE: In current version, LAPACK is used)
# -framework Accelerate
FFLAGS  := -O3 -framework Accelerate \
           -fbacktrace -ffpe-summary=none \
           -fopenmp -m64 -march=native -flto \
           -I$(MOD) -J$(MOD)

##########################################################################################
# Do NOT edit below
##########################################################################################

exclude := iso_c_binding omp_lib
sources := $(notdir $(shell ls $(SRC)*.f03))

# Extract dependency from the 'use' sequence in source files
#   e.g. '  use constants, only i4b' -> 'constants'
# $1: %.f03
dependencies=$(filter-out $(exclude), \
               $(sort \
                 $(shell grep -E '^\s*use\s+' "$(addprefix $(SRC),$1)" | \
                         awk '{ print $$2 }' | \
                         sed 's/,//g' \
                  )\
                )\
              )

# Show dependencies
tree=$(foreach file, $(sources), \
	     $(shell ls ./) \
	    )
#tree=$(foreach file, $(sources), \
#	     $(shell echo "$(file) <- $(call dependencies, $(file))\n") \
#	    )

# Rule for object files; %.o: %.f03 and its dependencies
# $1: $(OBJ)%.o
# $2: $(SRC)%.f03
define object
$1: $2 \
    $(patsubst %,$(OBJ)%.o, \
      $(call dependencies, \
        $(patsubst %.o,%.f03, \
          $(notdir $1)\
         )\
       )\
     )
	@mkdir -p $(OBJ)
	@mkdir -p $(MOD)
	$(FC) $(FFLAGS) -c $2 -o $$@
endef

# Main
$(BIN)$(PROGRAM): $(patsubst %.f03,$(OBJ)%.o,$(sources))
	@mkdir -p $(BIN)
	$(FC) $(FFLAGS) $^ -o $@

# Define object rules for all sources
$(foreach file, $(sources), \
  $(eval \
    $(call object, \
      $(patsubst %.f03, $(OBJ)%.o, $(file)),\
      $(patsubst %.f03, $(SRC)%.f03, $(file))\
     )\
   )\
 )

$(MOD)%.mod: $(SRC)%.f03 $(OBJ)%.o
	@:

.PHONY: clean install debug all tree
all: clean $(BIN)$(PROGRAM)
	@:

clean:
	- rm -rf $(BIN) $(OBJ) $(MOD)

install: $(BIN)$(PROGRAM)
	install -s $(BIN)$(PROGRAM) $(DEST)

# Show dependencies (not doing anything else)
tree:
	@:
	@$(foreach file,$(sources), $(warning $(shell echo $(file) | awk '{ printf "%20s", $$1 }') <- $(call dependencies,$(file))))

debug:
	@echo $(sources)
	@echo $(patsubst %.f03,$(OBJ)%.o,$(sources))
	@echo $(patsubst %,$(OBJ)%.o,$(call dependencies,bz.f03))


