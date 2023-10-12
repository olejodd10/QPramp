CC := gcc
FLAGS := -Ofast -g

BUILD_DIR := build
EXECUTABLE := main
LIB := qpramp
ARCHIVE := $(BUILD_DIR)/lib$(LIB).a
SOURCES := qp_ramp.c csv.c timing.c vector.c matrix.c lti.c iterable_set.c
OBJECTS := $(patsubst %.c, $(BUILD_DIR)/%.o, $(SOURCES)) 

MEX_DIR := mex
MEX_FUNCTIONS := qp_ramp_solve_mpc qp_ramp_solve
MEX_FILES := $(patsubst %, $(MEX_DIR)/%.mexa64, $(MEX_FUNCTIONS)) 
MATLAB_PATH := /home/ole/programs/matlab/bin/matlab

.PHONY: all run clean lib mex
all: run

run: $(EXECUTABLE)
	./$<

lib: $(ARCHIVE)

$(BUILD_DIR):
	mkdir $@

$(BUILD_DIR)/%.o: %.c | $(BUILD_DIR)
	$(CC) $(FLAGS) -c $< -o $@

$(ARCHIVE): $(OBJECTS)
	ar rcs $@ $(OBJECTS)

$(EXECUTABLE): $(EXECUTABLE).c $(ARCHIVE)
	$(CC) $(FLAGS) $< -o $@ -L. -l:$(ARCHIVE)

mex: $(MEX_FILES)

$(MEX_DIR)/%.mexa64: $(MEX_DIR)/%.c $(ARCHIVE)
	$(MATLAB_PATH) -nodesktop -nosplash -r "mex -L. -l:$(ARCHIVE) -outdir $(MEX_DIR) $<; exit;"

clean:
	rm -rf $(BUILD_DIR)
