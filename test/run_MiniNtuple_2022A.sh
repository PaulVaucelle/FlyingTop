#!/bin/bash

# Compile executable
printf "\n>>> Compiling MiniNtuple_2022...\n"
COMPILER=$(root-config --cxx)
FLAGS=$(root-config --cflags --libs)
$COMPILER -g -O3 -o MiniNtuple_RunSignal2022  MiniNtuple_RunSignal2022A.C MiniNtuple.C $FLAGS

printf "\n>>> Running MiniNtuple_S2022A\n"
./MiniNtuple_RunSignal2022
