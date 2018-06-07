#!/bin/bash

mkdir ../../Data/simout070718_9
mkdir ../../Data/simout070718_10
mkdir ../../Data/simout070718_11
mkdir ../../Data/simout070718_12

./L1-cellpop-lite.r batch 17 ../../Data/simout070718_9/out &
echo $! > ../../Data/simout070718_9/pid.txt
echo "./L1-cellpop-lite.r batch 17 ../../Data/simout070718_9 &" >> ../../Data/simout070718_9/pid.txt

./L1-cellpop-lite.r batch 17 ../../Data/simout070718_10/out &
echo $! > ../../Data/simout070718_10/pid.txt
echo "./L1-cellpop-lite.r batch 17 ../../Data/simout070718_10 &" >> ../../Data/simout070718_10/pid.txt

./L1-cellpop-lite.r batch 17 ../../Data/simout070718_11/out &
echo $! > ../../Data/simout070718_11/pid.txt
echo "./L1-cellpop-lite.r batch 17 ../../Data/simout070718_11 &" >> ../../Data/simout070718_11/pid.txt

./L1-cellpop-lite.r batch 17 ../../Data/simout070718_12/out &
echo $! > ../../Data/simout070718_12/pid.txt
echo "./L1-cellpop-lite.r batch 17 ../../Data/simout070718_12 &" >> ../../Data/simout070718_12/pid.txt
