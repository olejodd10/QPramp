CONFIG := config.Makefile
include $(CONFIG)

INPUT_CSVS := $(patsubst %, $(INPUT_DIR_RAW)/%.csv, $(INPUT_VECTORS) $(INPUT_MATRICES))
LIB_PATHS := /home/ole/repos/qp-ramp/qpOASES-releases-3.2.1/bin
LIBS := -l:libqpOASES_wrapper.a -lstdc++ -lm
LIB_FLAGS := $(addprefix -L, $(LIB_PATHS)) $(LIBS) 
INCLUDE_PATHS := /home/ole/repos/qp-ramp/qpOASES-releases-3.2.1/interfaces/c
INCLUDE_FLAGS := $(addprefix -I, $(INCLUDE_PATHS))
DEFINE_FLAGS := $(foreach val, $(DEFINES), -D$(val)=$($(val)))

CC := gcc
FLAGS := -Ofast $(DEFINE_FLAGS) $(INCLUDE_FLAGS) $(LIB_FLAGS)
BUILD_DIR := build
EXECUTABLE := main
SOURCES := algs.c csv.c timing.c vector.c matrix.c lti.c iterable_set.c $(EXECUTABLE).c
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
	$(CC) -c $< -o $@ $(FLAGS) 

$(EXECUTABLE): $(OBJECTS) $(CONFIG)
	$(CC) $(OBJECTS) -o $@ $(FLAGS) 

$(SAVE_EXECUTABLE) : $(SAVE_OBJECTS) $(CONFIG)
	$(CC) $(SAVE_OBJECTS) -o $@ $(FLAGS) 

$(SAVE_CSVS): $(SAVE_EXECUTABLE) $(INPUT_CSVS)
	./$<

runsave: $(SAVE_CSVS)

clean:
	rm -rf $(BUILD_DIR)
	rm $(EXECUTABLE)

cleansave:
	rm $(SAVE_CSVS)
	rm $(SAVE_EXECUTABLE)
