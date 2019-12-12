#!/bin/bash

echo "Compiling CUDA400 ..."
make -C CUDA400

echo "Compiling CUDA512 ..."
make -C CUDA512

echo "Compiling CUDA1024 ..."
make -C CUDA1024

echo "Compiling CUDA2048 ..."
make -C CUDA2048

echo "Compiling CUDA4096 ..."
make -C CUDA4096

echo "Compiling CUDA8192 ..."
make -C CUDA8192

#echo "Compiling CUDA16384 ..."
#make -C CUDA16384

#echo "Compiling CUDA32768 ..."
#make -C CUDA32768

#echo "Compiling CUDA65536 ..."
#make -C CUDA65536

#echo "Compiling CUDA131072 ..."
#make -C CUDA131072


