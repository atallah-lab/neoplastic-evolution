#!/bin/bash

./L1-cellpop-lite-v22.r batch ../../sims080818/run 50 0 0.1 0.001 &
echo $! > ../../sims080818/pid.txt
./L1-cellpop-lite-v22.r batch ../../sims080818/run 100 0 0.1 0.001 &
echo $! >> ../../sims080818/pid.txt
./L1-cellpop-lite-v22.r batch ../../sims080818/run 500 0 0.1 0.001 &
echo $! >> ../../sims080818/pid.txt
./L1-cellpop-lite-v22.r batch ../../sims080818/run 1000 0 0.1 0.001 &
echo $! >> ../../sims080818/pid.txt
./L1-cellpop-lite-v22.r batch ../../sims080818/run 2000 0 0.1 0.001 &
echo $! >> ../../sims080818/pid.txt
