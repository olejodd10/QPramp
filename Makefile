CC = gcc
BUILDDIR = build
SOURCES = $(wildcard *.c)
OBJECTS = $(patsubst %.c, $(BUILDDIR)/%.o, $(SOURCES)) 
EXECUTABLE = main

.PHONY: all run clean
all: run

run: $(EXECUTABLE)
	./$<

$(BUILDDIR):
	mkdir $@

$(BUILDDIR)/%.o: %.c | $(BUILDDIR)
	$(CC) -c $< -o $@

$(EXECUTABLE): $(OBJECTS)
	$(CC) $^ -o $(EXECUTABLE)

clean:
	rm -rf $(BUILDDIR)
	rm $(EXECUTABLE)
