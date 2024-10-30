#!/bin/bash

# Compile executable
printf "\n>>> Compiling MiniNtuple_2023...\n"
COMPILER=$(root-config --cxx)
FLAGS=$(root-config --cflags --libs)
$COMPILER -g -O3 -o MiniNtuple_2023  MiniNtuple_RunData2023.C MiniNtuple.C $FLAGS

printf "\n>>> Running MiniNtuple_2023\n"
./MiniNtuple_2023
