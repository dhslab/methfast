# Makefile for methfast program

# Compiler
CC = gcc

# Flags for compilation
CFLAGS = -O2 -Wall

# Source files
SOURCES = methfast.c cgranges.c

# Object files
OBJECTS = $(SOURCES:.c=.o)

# Executable name
EXEC = methfast

# Installation directory
INSTALL_DIR = /usr/bin

# Default target: build the program
all: $(EXEC)

# Build the executable
$(EXEC): $(OBJECTS)
	$(CC) $(CFLAGS) -o $(EXEC) $(OBJECTS) -lz -lm

# Compile .c files to .o files
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Clean up object files and executable
clean:
	rm -f $(OBJECTS) $(EXEC)

# Install the executable to /usr/bin
install: $(EXEC)
	install -m 755 $(EXEC) $(INSTALL_DIR)

# Uninstall the program
uninstall:
	rm -f $(INSTALL_DIR)/$(EXEC)

.PHONY: all clean install uninstall
