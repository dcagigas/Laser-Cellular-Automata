#!/bin/bash

#rm time_results.txt

today=`date +%Y-%m-%d.%H:%M:%S`


echo "CUDA400 ..."
cd CUDA400
echo " " >> time_results400.txt
echo "CUDA Laser Simulation times: "$today >> time_results400.txt
echo "----------------------------" >> time_results400.txt
echo " " >> time_results400.txt
X1=`(time ./laserCAv6_400) 2>&1 | grep real | cut -c5-20`
pid="$!"
wait $pid 
echo '400:' $X1 >> time_results400.txt
cd ..


echo "CUDA512 ..."
cd CUDA512
echo " " >> time_results512.txt
echo "CUDA Laser Simulation times: "$today >> time_results512.txt
echo "----------------------------" >> time_results512.txt
echo " " >> time_results512.txt
X2=`(time ./laserCAv6_512) 2>&1 | grep real | cut -c5-20`
pid="$!"
wait $pid 
echo '512:' $X2 >> time_results512.txt
cd ..


echo "CUDA1024 ..."
cd CUDA1024
echo " " >> time_results1024.txt
echo "CUDA Laser Simulation times: "$today >> time_results1024.txt
echo "----------------------------" >> time_results1024.txt
echo " " >> time_results1024.txt
X3=`(time ./laserCAv6_1024) 2>&1 | grep real | cut -c5-20`
pid="$!"
wait $pid 
echo '1024:' $X3 >> time_results1024.txt
cd ..


echo "CUDA2048 ..."
cd CUDA2048
echo " " >> time_results2048.txt
echo "CUDA Laser Simulation times: "$today >> time_results2048.txt
echo "----------------------------" >> time_results2048.txt
echo " " >> time_results2048.txt
X4=`(time ./laserCAv6_2048) 2>&1 | grep real | cut -c5-20`
pid="$!"
wait $pid 
echo '2048:' $X4 >> time_results2048.txt
cd ..


echo "CUDA4096 ..."
cd CUDA4096
echo " " >> time_results4096.txt
echo "CUDA Laser Simulation times: "$today >> time_results4096.txt
echo "----------------------------" >> time_results4096.txt
echo " " >> time_results4096.txt
X5=`(time ./laserCAv6_4096) 2>&1 | grep real | cut -c5-20`
pid="$!"
wait $pid 
echo '4096:' $X5 >> time_results4096.txt
cd ..


#echo "CUDA8192 ..."
#cd CUDA8192
#echo " " >> time_results8192.txt
#echo "CUDA Laser Simulation times: "$today >> time_results8192.txt
#echo "----------------------------" >> time_results8192.txt
#echo " " >> time_results8192.txt
#X6=`(time ./laserCAv6_8192) 2>&1 | grep real | cut -c5-20`
#pid="$!"
#wait $pid 
#echo '8192:' $X6 >> time_results8192.txt
#cd ..


