#!/bin/bash

mkdir ../../Data/simout050718_1
mkdir ../../Data/simout050718_2
mkdir ../../Data/simout050718_3
mkdir ../../Data/simout050718_4
#mkdir ../../Data/simout050718_5

./L1-cellpop.r batch 5 ../../Data/simout050718_1/out &
echo $! > ../../Data/simout050718_1/pid.txt
echo "./L1-cellpop.r batch 5 ../../Data/simout050718_1 &" >> ../../Data/simout050718_1/pid.txt

./L1-cellpop.r batch 5 ../../Data/simout050718_2/out &
echo $! > ../../Data/simout050718_2/pid.txt
echo "./L1-cellpop.r batch 5 ../../Data/simout050718_2 &" >> ../../Data/simout050718_2/pid.txt

./L1-cellpop.r batch 5 ../../Data/simout050718_3/out &
echo $! > ../../Data/simout050718_3/pid.txt
echo "./L1-cellpop.r batch 5 ../../Data/simout050718_3 &" >> ../../Data/simout050718_3/pid.txt

./L1-cellpop.r batch 5 ../../Data/simout050718_4/out &
echo $! > ../../Data/simout050718_4/pid.txt
echo "./L1-cellpop.r batch 5 ../../Data/simout050718_4 &" >> ../../Data/simout050718_4/pid.txt

#./L1-cellpop.r batch 5 ../../Data/simout050718_5 &
#echo $! > ../../Data/simout050718_5
#echo "./L1-cellpop.r batch 5 ../../Data/simout050718_5 &" >> ../../Data/simout050718_5
