#!/bin/bash

# Compile executable
printf "\n>>> Compiling MiniDATAMCNtuple_2023...\n"
COMPILER=$(root-config --cxx)
FLAGS=$(root-config --cflags --libs)
$COMPILER -g -O3 -o MiniDATAMCNtuple_RunData2023B  MiniDATAMCNtuple_RunData2023B.C MiniDATAMCNtuple.C $FLAGS

printf "\n>>> Running MiniNtuple_RunData2023B\n"
./MiniDATAMCNtuple_RunData2023B
