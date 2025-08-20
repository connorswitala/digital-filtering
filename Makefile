# ==== Master Makefile ====

# Subdirectories
CPP_DIR      := digital-filtering-c++
FORTRAN_DIR  := digital-filtering-fortran

# Default target: build both
all: cpp fortran

# Build C++
cpp:
	$(MAKE) -C $(CPP_DIR)

# Build Fortran
fortran:
	$(MAKE) -C $(FORTRAN_DIR)

# Clean only C++ outputs
clean-cpp:
	$(MAKE) -C $(CPP_DIR) clean

# Clean only Fortran outputs
clean-fortran:
	$(MAKE) -C $(FORTRAN_DIR) clean

# Clean everything
clean: clean-cpp clean-fortran

.PHONY: all cpp fortran clean clean-cpp clean-fortran
