#!/bin/bash

mkdir ../../Data/simout053018_9
mkdir ../../Data/simout053018_10
mkdir ../../Data/simout053018_11
mkdir ../../Data/simout053018_12

./L1-cellpop-lite.r batch 17 ../../Data/simout053018_9/out &
echo $! > ../../Data/simout053018_9/pid.txt
echo "./L1-cellpop-lite.r batch 17 ../../Data/simout053018_9 &" >> ../../Data/simout053018_9/pid.txt

./L1-cellpop-lite.r batch 17 ../../Data/simout053018_10/out &
echo $! > ../../Data/simout053018_10/pid.txt
echo "./L1-cellpop-lite.r batch 17 ../../Data/simout053018_10 &" >> ../../Data/simout053018_10/pid.txt

./L1-cellpop-lite.r batch 17 ../../Data/simout053018_11/out &
echo $! > ../../Data/simout053018_11/pid.txt
echo "./L1-cellpop-lite.r batch 17 ../../Data/simout053018_11 &" >> ../../Data/simout053018_11/pid.txt

./L1-cellpop-lite.r batch 17 ../../Data/simout053018_12/out &
echo $! > ../../Data/simout053018_12/pid.txt
echo "./L1-cellpop-lite.r batch 17 ../../Data/simout053018_12 &" >> ../../Data/simout053018_12/pid.txt
