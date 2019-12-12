#!/bin/bash

today=`date +%Y-%m-%d.%H:%M:%S`


echo "OPENMP400 ..."
cd OPENMP400
echo " " >> time_results400.txt
echo "CUDA Laser Simulation times: "$today >> time_results400.txt
echo "----------------------------" >> time_results400.txt
echo " " >> time_results400.txt
X1=`(time ./laseracv10_omp_basic_V0) 2>&1 | grep real | cut -c5-20`
pid="$!"
wait $pid 
echo '400:' $X1 >> time_results400.txt
cd ..


echo "OPENMP512 ..."
cd OPENMP512
echo " " >> time_results512.txt
echo "CUDA Laser Simulation times: "$today >> time_results512.txt
echo "----------------------------" >> time_results512.txt
echo " " >> time_results512.txt
X2=`(time ./laseracv10_omp_basic_V0) 2>&1 | grep real | cut -c5-20`
pid="$!"
wait $pid 
echo '512:' $X2 >> time_results512.txt
cd ..


echo "OPENMP1024 ..."
cd OPENMP1024
echo " " >> time_results1024.txt
echo "CUDA Laser Simulation times: "$today >> time_results1024.txt
echo "----------------------------" >> time_results1024.txt
echo " " >> time_results1024.txt
X3=`(time ./laseracv10_omp_basic_V0) 2>&1 | grep real | cut -c5-20`
pid="$!"
wait $pid 
echo '1024:' $X3 >> time_results1024.txt
cd ..


echo "OPENMP2048 ..."
cd OPENMP2048
echo " " >> time_results2048.txt
echo "CUDA Laser Simulation times: "$today >> time_results2048.txt
echo "----------------------------" >> time_results2048.txt
echo " " >> time_results2048.txt
X4=`(time ./laseracv10_omp_basic_V0) 2>&1 | grep real | cut -c5-20`
pid="$!"
wait $pid 
echo '2048:' $X4 >> time_results2048.txt
cd ..


echo "OPENMP4096 ..."
cd OPENMP4096
echo " " >> time_results4096.txt
echo "CUDA Laser Simulation times: "$today >> time_results4096.txt
echo "----------------------------" >> time_results4096.txt
echo " " >> time_results4096.txt
X5=`(time ./laseracv10_omp_basic_V0) 2>&1 | grep real | cut -c5-20`
pid="$!"
wait $pid 
echo '4096:' $X5 >> time_results4096.txt
cd ..


