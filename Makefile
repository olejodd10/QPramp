LIB_PATHS := /home/ole/repos/qp-ramp/BLAS/OpenBLAS/lib /home/ole/repos/qp-ramp/osqp-v0.6.3-Linux/lib
LIBS := -l:libosqp.a -l:libopenblas.a -lm 
LIB_FLAGS := $(addprefix -L, $(LIB_PATHS)) $(LIBS) 
INCLUDE_PATHS := /home/ole/repos/qp-ramp/BLAS/OpenBLAS/include /home/ole/repos/qp-ramp/osqp-v0.6.3-Linux/include
INCLUDE_FLAGS := $(addprefix -I, $(INCLUDE_PATHS))

CC := gcc
FLAGS := -Ofast -g $(INCLUDE_FLAGS) $(LIB_FLAGS)

BUILD_DIR := build
EXECUTABLE := main
LIB := qpramp
ARCHIVE := $(BUILD_DIR)/lib$(LIB).a
SOURCES := qp_ramp.c csv.c timing.c vector.c matrix.c lti.c iterable_set.c osqp_init.c
OBJECTS := $(patsubst %.c, $(BUILD_DIR)/%.o, $(SOURCES)) 

.PHONY: all run clean lib
all: run

run: $(EXECUTABLE)
	./$<

lib: $(ARCHIVE)

$(BUILD_DIR):
	mkdir $@

$(BUILD_DIR)/%.o: %.c | $(BUILD_DIR)
	$(CC) -c $< -o $@ $(FLAGS) 

$(ARCHIVE): $(OBJECTS)
	ar rcs $@ $(OBJECTS)

$(EXECUTABLE): $(EXECUTABLE).c $(ARCHIVE)
	$(CC) $< -o $@ -L. -l:$(ARCHIVE) $(FLAGS) 

clean:
	rm -rf $(BUILD_DIR)
