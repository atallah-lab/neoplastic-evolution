#!/bin/bash

./L1-cellpop-lite-v22.r batch ../../sims062718/run 50 0 0 0 > ../../sims062718/run &
echo $! > ../../sims062718/pid.txt
./L1-cellpop-lite-v22.r batch ../../sims062718/run 100 0 0 0 > ../../sims062718/run &
echo $! >> ../../sims062718/pid.txt
#./L1-cellpop-lite-v22.r batch ../../sims062718/run 500 0 0 0 > ../../sims062718/run &
#echo $! >> ../../sims062718/pid.txt
#./L1-cellpop-lite-v22.r batch ../../sims062718/run 1000 0 0 0 > ../../sims062718/run &
#echo $! >> ../../sims062718/pid.txt
#./L1-cellpop-lite-v22.r batch ../../sims062718/run 2000 0 0 0 > ../../sims062718/run &
#echo $! >> ../../sims062718/pid.txt
