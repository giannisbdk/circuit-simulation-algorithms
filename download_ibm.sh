#!/bin/sh

if [ -d "cd ibm_netlists" ]; then
	cd ibm_netlists
else
	mkdir ibm_netlists
	cd ibm_netlists
fi

FILE=$1     

echo "\n"

case $FILE in
	"ibm1" )
		if [ -f "ibmpg1.spice" ]; then
		   echo "File ibmpg1.spice already exists!"
		else
			wget http://dropzone.tamu.edu/~pli/PGBench/ibmpg/ibmpg1.spice.bz2
			bzip2 -d ibmpg1.spice.bz2 
		fi
		if [ -f "ibmpg1.solution" ]; then
		   echo "File ibmpg1.solution already exists!"
		else
			wget http://dropzone.tamu.edu/~pli/PGBench/ibmpg/ibmpg1.solution.bz2
			bzip2 -d ibmpg1.solution.bz2
		fi
	;;
	"ibm2" )
		if [ -f "ibmpg2.spice" ]; then
		   echo "File ibmpg2.spice already exists!"
		else
			wget http://dropzone.tamu.edu/~pli/PGBench/ibmpg/ibmpg2.spice.bz2
			bzip2 -d ibmpg2.spice.bz2 
		fi
		if [ -f "ibmpg2.solution" ]; then
		   echo "File ibmpg2.solution already exists!"
		else
			wget http://dropzone.tamu.edu/~pli/PGBench/ibmpg/ibmpg2.solution.bz2
			bzip2 -d ibmpg2.solution.bz2
		fi
	;;
	"ibm3" )
		if [ -f "ibmpg3.spice" ]; then
		   echo "File ibmpg3.spice already exists!"
		else
			wget http://dropzone.tamu.edu/~pli/PGBench/ibmpg/ibmpg3.spice.bz2
			bzip2 -d ibmpg3.spice.bz2 
		fi
		if [ -f "ibmpg3.solution" ]; then
		   echo "File ibmpg3.solution already exists!"
		else
			wget http://dropzone.tamu.edu/~pli/PGBench/ibmpg/ibmpg3.solution.bz2
			bzip2 -d ibmpg3.solution.bz2
		fi
	;;
	"ibm4" )
		if [ -f "ibmpg4.spice" ]; then
		   echo "File ibmpg4.spice already exists!"
		else
			wget http://dropzone.tamu.edu/~pli/PGBench/ibmpg/ibmpg4.spice.bz2
			bzip2 -d ibmpg4.spice.bz2 
		fi
		if [ -f "ibmpg4.solution" ]; then
		   echo "File ibmpg4.solution already exists!"
		else
			wget http://dropzone.tamu.edu/~pli/PGBench/ibmpg/ibmpg4.solution.bz2
			bzip2 -d ibmpg4.solution.bz2
		fi 
	;;
	"ibm5" )
		if [ -f "ibmpg5.spice" ]; then
		   echo "File ibmpg5.spice already exists!"
		else
			wget http://dropzone.tamu.edu/~pli/PGBench/ibmpg/ibmpg5.spice.bz2
			bzip2 -d ibmpg5.spice.bz2 
		fi
		if [ -f "ibmpg5.solution" ]; then
		   echo "File ibmpg5.solution already exists!"
		else
			wget http://dropzone.tamu.edu/~pli/PGBench/ibmpg/ibmpg5.solution.bz2
			bzip2 -d ibmpg5.solution.bz2
		fi 
	;;
	"ibm6" )
		if [ -f "ibmpg6.spice" ]; then
		   echo "File ibmpg6.spice already exists!"
		else
			wget http://dropzone.tamu.edu/~pli/PGBench/ibmpg/ibmpg6.spice.bz2
			bzip2 -d ibmpg6.spice.bz2 
		fi
		if [ -f "ibmpg6.solution" ]; then
		   echo "File ibmpg6.solution already exists!"
		else
			wget http://dropzone.tamu.edu/~pli/PGBench/ibmpg/ibmpg6.solution.bz2
			bzip2 -d ibmpg6.solution.bz2
		fi 
	;;
	*)
		echo "sorry did't understand"
esac

cd .. 

echo "\n\nREMEMBER TO SET THE APPROPRIATE OPTIONS IN THE SPICE FILES\n\n"