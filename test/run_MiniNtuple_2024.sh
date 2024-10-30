#!/bin/bash

# Compile executable
printf "\n>>> Compiling MiniNtuple_2024...\n"
COMPILER=$(root-config --cxx)
FLAGS=$(root-config --cflags --libs)
$COMPILER -g -O3 -o MiniNtuple_2024  MiniNtuple_RunData2024.C MiniNtuple.C $FLAGS

printf "\n>>> Running MiniNtuple_2024\n"
./MiniNtuple_2024
