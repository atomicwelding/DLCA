#!/bin/bash
set -e

THREADS=8
FLAGS="-O3 -march=native -fopenmp -Accelerate"
TARGET="dlca"

echo "This software is provided by weld"
echo "See https://github.com/atomicwelding"


echo "Forcing to rebuild ..."
find . -name "*.o" -delete
find . -name "*.mod" -delete


echo "Starting to build ..."

gfortran -c params.f90  $FLAGS
gfortran -c types.f90   $FLAGS
gfortran -c subroutines.f90 $FLAGS
gfortran -c main.f90    $FLAGS

echo "Linking ..."
gfortran types.o params.o subroutines.o main.o $FLAGS -o $TARGET

export OMP_NUM_THREADS=$THREADS
echo "Executable is ready. Running on $THREADS thread(s) : ./$TARGET"
./$TARGET
