# Project directories
SRC_DIR      := src
BUILD_DIR    := build
INCLUDE_DIR  := include

# Compiler and flags
FC           := gfortran
FFLAGS       := -g -fcheck=all -fbacktrace -Wall -Wextra -J$(INCLUDE_DIR)
LDFLAGS      :=
EXEC         := $(BUILD_DIR)/maniac

# Source and object files
SRCS         := $(wildcard $(SRC_DIR)/*.f90)
OBJS         := $(patsubst $(SRC_DIR)/%.f90, $(BUILD_DIR)/%.o, $(SRCS))

# Default target
all: $(EXEC)

# Link executable
$(EXEC): $(OBJS) | $(BUILD_DIR)
	$(FC) $(LDFLAGS) -o $@ $(OBJS)
	@echo "Build complete: $@"

# Compile source files to object files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90 | $(BUILD_DIR) $(INCLUDE_DIR)
	$(FC) $(FFLAGS) -c -o $@ $<

# Ensure build and include directories exist
$(BUILD_DIR) $(INCLUDE_DIR):
	mkdir -p $@

# Generate dependencies using makedepf90
depend: $(SRCS)
	makedepf90 -b $(BUILD_DIR) $(SRCS) > dependencies.d

# Include generated dependencies
-include dependencies.d

# Clean targets
clean:
	rm -rf $(BUILD_DIR)/*.o $(BUILD_DIR)/*.mod $(EXEC)

distclean: clean
	rm -rf $(INCLUDE_DIR) dependencies.d

.PHONY: all clean distclean depend

