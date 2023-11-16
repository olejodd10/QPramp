SRC_DIR := src
BUILD_DIR := build

INCLUDE_PATHS := include
INCLUDE_FLAGS := $(addprefix -I, $(INCLUDE_PATHS))

LIB_PATHS := 
LIBS := -lm -lpthread
LIB_FLAGS := $(addprefix -L, $(LIB_PATHS)) $(LIBS)

ARCHIVE := libqpramp.a
SOURCES := qp_ramp.c vector.c matrix.c iterable_set.c thread_pool.c parallelism.c
OBJECTS := $(patsubst %.c, $(BUILD_DIR)/%.o, $(SOURCES)) 

EXAMPLE := main
EXAMPLE_SOURCES := main.c lti.c csv.c timing.c matrix.c  
EXAMPLE_OBJECTS := $(patsubst %.c, $(BUILD_DIR)/%.o, $(EXAMPLE_SOURCES)) 

MEX_DIR := mex
MEX_FUNCTIONS := qp_ramp_solve_mpc qp_ramp_solve
MEX_FILES := $(patsubst %, $(MEX_DIR)/%.mexa64, $(MEX_FUNCTIONS)) 
MATLAB_PATH := /home/ole/programs/matlab/bin/matlab

CC := gcc
FLAGS := -Ofast -g

.PHONY: all clean lib mex
all: $(EXAMPLE)

lib: $(ARCHIVE)

$(BUILD_DIR):
	mkdir $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c | $(BUILD_DIR)
	$(CC) -c $< -o $@ $(FLAGS) $(INCLUDE_FLAGS)

$(ARCHIVE): $(OBJECTS)
	ar rcs $@ $(OBJECTS)

$(EXAMPLE): $(EXAMPLE_OBJECTS) $(ARCHIVE)
	$(CC) $(EXAMPLE_OBJECTS) -o $@ $(FLAGS) $(INCLUDE_FLAGS) -L. -l:$(ARCHIVE) $(LIB_FLAGS)

mex: $(MEX_FILES)

$(MEX_DIR)/%.mexa64: $(MEX_DIR)/%.c $(ARCHIVE)
	$(MATLAB_PATH) -nodesktop -nosplash -r "mex $(INCLUDE_FLAGS) -L. -l:$(ARCHIVE) $(LIB_FLAGS) -outdir $(MEX_DIR) $<; exit;"

clean:
	rm -rf $(BUILD_DIR)
	rm $(ARCHIVE)
