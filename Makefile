# Compiler
CC = gcc
# Compiler flags
CFLAGS = -g -Wall
# Compile with -O3 optimization
OPTFLAGS = -O3
# Directories
SRC = src
OBJ = obj

SOURCES = $(wildcard $(SRC)/*.c)
OBJECTS = $(patsubst $(SRC)/%.c,$(OBJ)/%.o,$(SOURCES))

main: $(OBJECTS)
	$(CC) $(CFLAGS) $^ -o $@

$(OBJ)/%.o: $(SRC)/%.c
	$(CC) $(CFLAGS) -I$(SRC) -c $< -o $@

.PHONY: clean
# Clean only the executable
clean:
	$(RM) main

# Cleans the executable and the object files
.PHONY: cleanall
cleanall: clean
	$(RM) $(OBJECTS)