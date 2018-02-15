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
IBM_NLS = ibm_netlists
PLT = plots

PLT_SRC = plot.py
CSX_LIB = cx_sparse/Lib/libcxsparse.a
CSX_SRC = cx_sparse/Source

SOURCES = $(wildcard $(SRC)/*.c)
OBJECTS = $(patsubst $(SRC)/%.c,$(OBJ)/%.o,$(SOURCES))

main: $(OBJECTS)
	$(CC) $(CFLAGS) $^ $(GSLFLAGS) -o $@ $(DBGFLAGS) $(CSX_LIB)

$(OBJ)/%.o: $(SRC)/%.c
	$(CC) $(CFLAGS) -I$(SRC) -c $< $(GSLFLAGS) -o $@ $(DBGFLAGS)

debug: DBGFLAGS = -DDEBUGL -DDEBUGH
debug: main

# Various targets with different input files
run_lu: main
	./main $(NLS)/lu_netlist.txt

run_lu_iter: main
	./main $(NLS)/lu_iter_netlist.txt

run_lu_sparse: main
	./main $(NLS)/lu_sparse_netlist.txt

run_lu_sparse_iter: main
	./main $(NLS)/lu_sparse_iter_netlist.txt

run_lu_tran_tr: main
	./main $(NLS)/lu_tran_tr_netlist.txt

run_lu_tran_tr_iter: main
	./main $(NLS)/lu_tran_tr_iter_netlist.txt

run_lu_tran_tr_sparse: main
	./main $(NLS)/lu_tran_tr_sparse_netlist.txt

run_lu_tran_tr_iter_sparse: main
	./main $(NLS)/lu_tran_tr_iter_sparse_netlist.txt

run_lu_tran_be: main
	./main $(NLS)/lu_tran_be_netlist.txt

run_lu_tran_be_iter: main
	./main $(NLS)/lu_tran_be_iter_netlist.txt

run_lu_tran_be_sparse: main
	./main $(NLS)/lu_tran_be_sparse_netlist.txt

run_lu_tran_be_iter_sparse: main
	./main $(NLS)/lu_tran_be_iter_sparse_netlist.txt

run_chol: main
	./main $(NLS)/cholesky_netlist.txt

run_chol_iter: main
	./main $(NLS)/cholesky_iter_netlist.txt

run_chol_sparse: main
	./main $(NLS)/cholesky_sparse_netlist.txt

run_chol_sparse_iter: main
	./main $(NLS)/cholesky_sparse_iter_netlist.txt

run_chol_tran_tr: main
	./main $(NLS)/cholesky_tran_tr_netlist.txt

run_chol_tran_tr_iter: main
	./main $(NLS)/cholesky_tran_tr_iter_netlist.txt

run_chol_tran_be: main
	./main $(NLS)/cholesky_tran_be_netlist.txt

run_chol_tran_be_iter: main
	./main $(NLS)/cholesky_tran_be_iter_netlist.txt

run_ac: main
	./main $(NLS)/ac_netlist.txt

run_ac_iter: main
	./main $(NLS)/ac_iter_netlist.txt

run_ac_sparse: main
	./main $(NLS)/ac_sparse_netlist.txt

run_ac_iter_sparse: main
	./main $(NLS)/ac_iter_sparse_netlist.txt

run_all: main
	./main $(NLS)/all_analyses_netlist.txt

# Netlist download targets
run_ibm1: main
	./main $(IBM_NLS)/ibmpg1.spice

run_ibm2: main
	./main $(IBM_NLS)/ibmpg2.spice

run_ibm3: main
	./main $(IBM_NLS)/ibmpg3.spice

run_ibm4: main
	./main $(IBM_NLS)/ibmpg4.spice

run_ibm5: main
	./main $(IBM_NLS)/ibmpg5.spice

run_ibm6: main
	./main $(IBM_NLS)/ibmpg6.spice

download_ibm1:
	./download_ibm.sh ibm1

download_ibm2:
	./download_ibm.sh ibm2

download_ibm3:
	./download_ibm.sh ibm3

download_ibm4:
	./download_ibm.sh ibm4

download_ibm5:
	./download_ibm.sh ibm5

download_ibm6:
	./download_ibm.sh ibm6	

# Plots the output files
plot:
	@echo "Executing $(PLT_SRC)..."
	@python2.7 $(PLT_SRC)

# Installs dependencies
deps:
	@echo "Installing dependencies..."
	sudo apt install libgsl0-dev
	sudo apt install python-matplotlib

.PHONY: clean
# Clean the output files
clean:
	$(RM) dc_*.txt tr_*.txt ac_*.txt

.PHONY: clean_ibm
# Clean the output files
clean_ibm:
	rm -rf ibm_netlists

.PHONY: clean_plots
# Clean the plots inside $(PLT) directory
clean_plots:
	$(RM) $(PLT)/*.png

# Cleans the executable, the output files and the object files
.PHONY: cleanall
cleanall: clean
	$(RM) main $(OBJECTS) -rf $(PLT)
