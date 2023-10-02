SIMULATION_TIMESTEPS := 100
EPS := 2.2204e-8

INPUT_DIR_RAW := ../examples/example3
INPUT_DIR := \"$(INPUT_DIR_RAW)\"
OUTPUT_DIR := \"$(INPUT_DIR_RAW)/out\"
REFERENCE_DIR := \"$(INPUT_DIR_RAW)/reference\"

SAVE_FORMAT := \"%.4le\"

DEFINES := SIMULATION_TIMESTEPS EPS INPUT_DIR OUTPUT_DIR SAVE_FORMAT REFERENCE_DIR
DEFINE_FLAGS := $(foreach val, $(DEFINES), -D$(val)=$($(val)))

CC := gcc
FLAGS := -Ofast $(DEFINE_FLAGS)

BUILD_DIR := build
EXECUTABLE := main
LIB := qpramp
ARCHIVE := $(BUILD_DIR)/lib$(LIB).a
SOURCES := qp_ramp.c csv.c timing.c vector.c matrix.c lti.c iterable_set.c
OBJECTS := $(patsubst %.c, $(BUILD_DIR)/%.o, $(SOURCES)) 

.PHONY: all run clean lib
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

clean:
	rm -rf $(BUILD_DIR)
