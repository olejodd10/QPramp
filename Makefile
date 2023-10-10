CC := gcc
FLAGS := -Ofast -g

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
