#!/bin/bash

echo "Removing compilated files OPENMP400 ..."
make -C OPENMP400 clean

echo "Removing compilated files OPENMP512 ..."
make -C OPENMP512 clean

echo "Removing compilated files OPENMP1024 ..."
make -C OPENMP1024 clean

echo "Removing compilated files OPENMP2048 ..."
make -C OPENMP2048 clean

echo "Removing compilated files OPENMP4096 ..."
make -C OPENMP4096 clean




