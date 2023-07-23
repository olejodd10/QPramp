CONFIG := config.Makefile
include $(CONFIG)

INPUT_CSVS := $(patsubst %, $(INPUT_DIR_RAW)/%.csv, $(INPUT_VECTORS) $(INPUT_MATRICES))
DEFINE_FLAGS := $(foreach val, $(DEFINES), -D$(val)=$($(val)))

CC := gcc
FLAGS := -Ofast $(DEFINE_FLAGS)
BUILD_DIR := build
SOURCES := $(wildcard *.c)
OBJECTS := $(patsubst %.c, $(BUILD_DIR)/%.o, $(SOURCES)) 
EXECUTABLE := main

.PHONY: all run clean
all: run

run: $(EXECUTABLE)
	./$<

$(BUILD_DIR):
	mkdir $@

$(BUILD_DIR)/%.o: %.c $(INPUT_CSVS) $(CONFIG) | $(BUILD_DIR)
	$(CC) $(FLAGS) -c $< -o $@

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(FLAGS) $^ -o $(EXECUTABLE)

clean:
	rm -rf $(BUILD_DIR)
	rm $(EXECUTABLE)
