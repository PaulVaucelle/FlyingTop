#!/bin/bash

# Compile executable
printf "\n>>> Compiling MiniNtuple_2022...\n"
COMPILER=$(root-config --cxx)
FLAGS=$(root-config --cflags --libs)
$COMPILER -g -O3 -o MiniNtuple_S2022B  MiniNtuple_RunSignal2022B.C MiniNtuple.C $FLAGS

printf "\n>>> Running MiniNtuple_S2022B\n"
./MiniNtuple_S2022B
