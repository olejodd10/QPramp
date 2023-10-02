CONFIG := config.Makefile
include $(CONFIG)

INPUT_CSVS := $(patsubst %, $(INPUT_DIR_RAW)/%.csv, $(INPUT_VECTORS) $(INPUT_MATRICES))
DEFINE_FLAGS := $(foreach val, $(DEFINES), -D$(val)=$($(val)))

CC := gcc
FLAGS := -Ofast $(DEFINE_FLAGS)
BUILD_DIR := build
EXECUTABLE := main
SOURCES := qp_ramp.c csv.c timing.c vector.c matrix.c lti.c iterable_set.c $(EXECUTABLE).c
OBJECTS := $(patsubst %.c, $(BUILD_DIR)/%.o, $(SOURCES)) 

SAVE_VECTORS := neg_w 
SAVE_MATRICES := neg_g_invh_gt neg_s neg_invh_f neg_g_invh  
SAVE_CSVS := $(patsubst %, $(INPUT_DIR_RAW)/%.csv, $(SAVE_VECTORS) $(SAVE_MATRICES))
SAVE_EXECUTABLE := save
SAVE_SOURCES := csv.c timing.c vector.c matrix.c $(SAVE_EXECUTABLE).c
SAVE_OBJECTS := $(patsubst %.c, $(BUILD_DIR)/%.o, $(SAVE_SOURCES)) 

.PHONY: all run runsave clean cleansave
all: run

run: $(EXECUTABLE)
	./$<

$(BUILD_DIR):
	mkdir $@

$(BUILD_DIR)/%.o: %.c $(CONFIG) | $(BUILD_DIR)
	$(CC) $(FLAGS) -c $< -o $@

$(EXECUTABLE): $(OBJECTS) $(CONFIG)
	$(CC) $(FLAGS) $(OBJECTS) -o $@

$(SAVE_EXECUTABLE) : $(SAVE_OBJECTS) $(CONFIG)
	$(CC) $(FLAGS) $(SAVE_OBJECTS) -o $@

$(SAVE_CSVS): $(SAVE_EXECUTABLE) $(INPUT_CSVS)
	./$<

runsave: $(SAVE_CSVS)

clean:
	rm -rf $(BUILD_DIR)
	rm $(EXECUTABLE)

cleansave:
	rm $(SAVE_CSVS)
	rm $(SAVE_EXECUTABLE)
