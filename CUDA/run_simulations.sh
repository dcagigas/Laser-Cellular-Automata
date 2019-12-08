#!/bin/bash

#rm time_results.txt

today=`date +%Y-%m-%d.%H:%M:%S`
echo " " >> time_results.txt
echo "CUDA Laser Simulation times: "$today >> time_results.txt
echo "----------------------------" >> time_results.txt
echo " " >> time_results.txt

X1=`(time ./laserCAv6_400 input/input400.dat) 2>&1 | grep real | cut -c5-20`
pid="$!"
wait $pid 
echo '400:' $X1 >> time_results.txt


X2=`(time ./laserCAv6_512 input/input512.dat) 2>&1 | grep real | cut -c5-20`
pid="$!"
wait $pid 
echo '512:' $X2 >> time_results.txt


X3=`(time ./laserCAv6_1024 input/input1024.dat) 2>&1 | grep real | cut -c5-20`
pid="$!"
wait $pid 
echo '1024:' $X3 >> time_results.txt


X4=`(time ./laserCAv6_2048 input/input2048.dat) 2>&1 | grep real | cut -c5-20`
pid="$!"
wait $pid
echo '2048:' $X4 >> time_results.txt


X5=`(time ./laserCAv6_4096 input/input4096.dat) 2>&1 | grep real | cut -c5-20`
pid="$!"
wait $pid
echo '4096:' $X5 >> time_results.txt


#X6=`(time ./laserCAv6_8192 input/input8192.dat) 2>&1 | grep real | cut -c5-20`
#pid="$!"
#wait $pid
#echo '8192:' $X6 >> time_results.txt

