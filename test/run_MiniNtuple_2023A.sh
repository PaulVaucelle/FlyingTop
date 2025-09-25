#!/bin/bash

# Compile executable
printf "\n>>> Compiling MiniNtuple_2023...\n"
COMPILER=$(root-config --cxx)
FLAGS=$(root-config --cflags --libs)
$COMPILER -g -O3 -o MiniNtuple_S2023A  MiniNtuple_RunSignal2023A.C MiniNtuple.C $FLAGS

printf "\n>>> Running MiniNtuple_S2023A\n"
./MiniNtuple_S2023A
