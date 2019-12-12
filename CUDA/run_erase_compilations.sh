#!/bin/bash

echo "Removing compilated files CUDA400 ..."
make -C CUDA400 clean

echo "Removing compilated files CUDA512 ..."
make -C CUDA512 clean

echo "Removing compilated files CUDA1024 ..."
make -C CUDA1024 clean

echo "Removing compilated files CUDA2048 ..."
make -C CUDA2048 clean

echo "Removing compilated files CUDA4096 ..."
make -C CUDA4096 clean

echo "Removing compilated files CUDA8192 ..."
make -C CUDA8192 clean

#echo "Removing compilated files CUDA16384 ..."
#make -C CUDA16384 clean

#echo "Removing compilated files CUDA32768 ..."
#make -C CUDA32768 clean

#echo "Removing compilated files CUDA65536 ..."
#make -C CUDA65536 clean

#echo "Removing compilated files CUDA131072 ..."
#make -C CUDA131072 clean



