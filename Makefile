SRC_DIR := src
BUILD_DIR := build

INCLUDE_PATHS := include
INCLUDE_FLAGS := $(addprefix -I, $(INCLUDE_PATHS))

LIB_PATHS := 
LIBS := -lm
LIB_FLAGS := $(addprefix -L, $(LIB_PATHS)) $(LIBS)

ARCHIVE := libqpramp.a
SOURCES := qp_ramp.c vector.c matrix.c iterable_set.c indexed_vectors.c
OBJECTS := $(patsubst %.c, $(BUILD_DIR)/%.o, $(SOURCES)) 

EXAMPLE := qp_ramp_example
EXAMPLE_SOURCES := main.c lti.c csv.c timing.c matrix.c  
EXAMPLE_OBJECTS := $(patsubst %.c, $(BUILD_DIR)/%.o, $(EXAMPLE_SOURCES)) 

CC := gcc
FLAGS := -Ofast -g -Wall

.PHONY: all clean lib
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

clean:
	rm -rf $(BUILD_DIR)
	rm $(ARCHIVE)
