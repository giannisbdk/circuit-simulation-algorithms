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

debug: DBGFLAGS = -DDEBUGL -DDEBUGH
debug: main

main: $(OBJECTS)
	$(CC) $(CFLAGS) $^ -o $@ $(DBGFLAGS)

$(OBJ)/%.o: $(SRC)/%.c
	$(CC) $(CFLAGS) -I$(SRC) -c $< -o $@ $(DBGFLAGS)

.PHONY: clean
# Clean only the executable
clean:
	$(RM) main

# Cleans the executable and the object files
.PHONY: cleanall
cleanall: clean
	$(RM) $(OBJECTS)