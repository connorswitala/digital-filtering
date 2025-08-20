# Compiler
FC = gfortran
FFLAGS = -O2 -Wall -fimplicit-none

# Directories
DF_DIR = df
TEST_DIR = test
BUILD_DIR = ../build/fortran
PROG_DIR = ../programs

# Sources
DF_SRC = $(DF_DIR)/df.f90
MAIN_SRC = $(TEST_DIR)/fortran-main.f90

# Objects
DF_OBJ = $(BUILD_DIR)/df.o
MAIN_OBJ = $(BUILD_DIR)/main.o

# Executable
EXE = $(PROG_DIR)/fortran-test.exe

# Default target
all: $(EXE)

# Build executable
$(EXE): $(DF_OBJ) $(MAIN_OBJ)
	@mkdir -p $(PROG_DIR)
	$(FC) $(DF_OBJ) $(MAIN_OBJ) -o $@

# Compile module
$(DF_OBJ): $(DF_SRC)
	@mkdir -p $(BUILD_DIR)
	$(FC) $(FFLAGS) -c $< -J $(BUILD_DIR) -o $@

# Compile main (depends on module)
$(MAIN_OBJ): $(MAIN_SRC) $(DF_OBJ)
	@mkdir -p $(BUILD_DIR)
	$(FC) $(FFLAGS) -c $< -I $(BUILD_DIR) -o $@

# Clean
clean:
	rm -f $(BUILD_DIR)/df.o $(BUILD_DIR)/main.o $(BUILD_DIR)/*.mod $(EXE)
