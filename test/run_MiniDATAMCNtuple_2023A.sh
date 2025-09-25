#!/bin/bash

# Compile executable
printf "\n>>> Compiling MiniDATAMCNtuple_2023...\n"
COMPILER=$(root-config --cxx)
FLAGS=$(root-config --cflags --libs)
$COMPILER -g -O3 -o MiniDATAMCNtuple_RunData2023A  MiniDATAMCNtuple_RunData2023A.C MiniDATAMCNtuple.C $FLAGS

printf "\n>>> Running MiniDATAMCNtuple_RunData2023A\n"
./MiniDATAMCNtuple_RunData2023A
