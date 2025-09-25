#!/bin/bash

# Compile executable
printf "\n>>> Compiling MiniDATAMCNtuple_2024...\n"
COMPILER=$(root-config --cxx)
FLAGS=$(root-config --cflags --libs)
$COMPILER -g -O3 -o MiniDATAMCNtuple_RunData2024  MiniDATAMCNtuple_RunData2024.C MiniDATAMCNtuple.C $FLAGS

printf "\n>>> Running MiniDATAMCNtuple_RunData2024\n"
./MiniDATAMCNtuple_RunData2024
