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
NLS = netlists

SOURCES = $(wildcard $(SRC)/*.c)
OBJECTS = $(patsubst $(SRC)/%.c,$(OBJ)/%.o,$(SOURCES))

main: $(OBJECTS)
	$(CC) $(CFLAGS) $^ $(GSLFLAGS) -o $@ $(DBGFLAGS)

$(OBJ)/%.o: $(SRC)/%.c
	$(CC) $(CFLAGS) -I$(SRC) -c $< $(GSLFLAGS) -o $@ $(DBGFLAGS)

debug: DBGFLAGS = -DDEBUGL -DDEBUGH
debug: main

# Various targets with different input files
run_lu: main
	./main $(NLS)/lu_netlist.txt

run_lu_iter: main
	./main $(NLS)/lu_iter_netlist.txt

run_lu_sparce: main
	./main $(NLS)/lu_sparse_netlist.txt

run_chol: main
	./main $(NLS)/cholesky_netlist.txt

run_chol_iter: main
	./main $(NLS)/cholesky_iter_netlist.txt

run_chol_sparce: main
	./main $(NLS)/cholesky_sparce_netlist.txt

.PHONY: clean
# Clean the output files
clean:
	$(RM) dc_*.txt

# Cleans the executable, the output files and the object files
.PHONY: cleanall
cleanall: clean
	$(RM) main $(OBJECTS)