# ==== Compiler and Flags ====
CXX       := g++
CXXFLAGS  := -std=c++17 -O2 

# ==== Directories ====
DF_DIR	:= df
TEST_DIR	:= test
BUILD_DIR	:= ../build/cpp
PROG_DIR	:= ../programs

# ==== Files ====
DF_SRC	:= $(DF_DIR)/df.cpp
DF_OBJ	:= $(BUILD_DIR)/df.o

TEST_SRC	:= $(TEST_DIR)/cpp-main.cpp
TEST_OBJ	:= $(BUILD_DIR)/main.o
TEST_TARGET	:= $(PROG_DIR)/cpp-test

# ==== Default target ====
all: $(TEST_TARGET)

# ==== Ensure build and programs dirs exist ====
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(PROG_DIR):
	mkdir -p $(PROG_DIR)

# ==== Object file compilation ====

$(TEST_OBJ): $(TEST_SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(DF_OBJ): $(DF_SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ==== Link targets ====	
$(TEST_TARGET): $(TEST_OBJ) $(DF_OBJ) | $(PROG_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@

# ==== Clean ====
clean:
	rm -f $(DF_OBJ) $(TEST_OBJ) $(TEST_TARGET)

.PHONY: all clean
