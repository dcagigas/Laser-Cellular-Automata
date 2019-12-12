#!/bin/bash

echo "Compiling OPENMP400 ..."
make -C OPENMP400

echo "Compiling OPENMP512 ..."
make -C OPENMP512

echo "Compiling OPENMP1024 ..."
make -C OPENMP1024

echo "Compiling OPENMP2048 ..."
make -C OPENMP2048

echo "Compiling OPENMP4096 ..."
make -C OPENMP4096

echo "Compiling OPENMP8192 ..."
make -C OPENMP8192

