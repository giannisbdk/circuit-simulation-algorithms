# Compiler
CC = gcc
# Compiler flags
CFLAGS = -g -Wall
# Compile with -O3 optimization
OPTFLAGS = -O3
# GSL config paths
GSLFLAGS = $(shell gsl-config --cflags --libs)
# Directories
SRC = src
OBJ = obj

SOURCES = $(wildcard $(SRC)/*.c)
OBJECTS = $(patsubst $(SRC)/%.c,$(OBJ)/%.o,$(SOURCES))

main: $(OBJECTS)
	$(CC) $(CFLAGS) $^ $(GSLFLAGS) -o $@ $(DBGFLAGS)

$(OBJ)/%.o: $(SRC)/%.c
	$(CC) $(CFLAGS) -I$(SRC) -c $< $(GSLFLAGS) -o $@ $(DBGFLAGS)

debug: DBGFLAGS = -DDEBUGL -DDEBUGH
debug: main

.PHONY: clean
# Clean only the executable
clean:
	$(RM) main

# Cleans the executable, the object files and the all the txt from analysis
.PHONY: cleanall
cleanall: clean
	$(RM) $(OBJECTS) dc_*.txt