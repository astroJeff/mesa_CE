

#################################################################

# STEP 1: set up the standard info

include $(MESA_DIR)/utils/makefile_header

#################################################################

# STEP 2: build the library

LIB = libbinary_CE.a

MOD_OTHER_DIR = $(MODULE_DIR)/other
MOD_DEFAULTS_DIR = $(MODULE_DIR)/defaults
OTHER_OBJS = other_binary_jdot.o

LIB_OBJS = binary_def.o binary_private_def.o\
   $(OTHER_OBJS) binary_jdot.o binary_mdot.o binary_tides.o binary_ctrls_io.o\
   binary_irradiation.o binary_history_specs.o binary_separation.o binary_history.o \
   binary_evolve.o binary_timestep.o\
   run_binary_support.o binary_lib.o

$(LIB) : $(LIB_OBJS)
	$(LIB_TOOL) $(LIB) $(LIB_OBJS)

BINARY_INCLUDES = -I$(MOD_DEFAULTS_DIR)
	
#################################################################

COMPILE = $(COMPILE_TO_DEPLOY) $(FCfree) $(BINARY_INCLUDES)
#COMPILE = $(COMPILE_TO_TEST) $(FCfree)

%.o: $(MOD_PUBLIC_DIR)/%.f
	$(COMPILE) $<

%.o: $(MOD_PRIVATE_DIR)/%.f
	$(COMPILE) $<

%.o: $(MOD_OTHER_DIR)/%.f
	$(COMPILE) $<

clean:
	-@rm -f *.o *.mod *.modmic *genmod.f90 *.a
