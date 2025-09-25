#!/bin/bash

# Compile executable
printf "\n>>> Compiling SecIntAna...\n"
COMPILER=$(root-config --cxx)
FLAGS=$(root-config --cflags --libs)
#$COMPILER -g -O3 -Wall -Wextra -Wpedantic -o SecIntAna SecIntAna_Run.C HistogramManager.C SecIntAna.C $FLAGS
$COMPILER -g -O3 -o SecIntAna SecIntAna_Run.C HistogramManager.C SecIntAna.C $FLAGS

printf "\n>>> Running SecIntAna\n"
./SecIntAna
