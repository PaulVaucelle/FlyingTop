#!/bin/bash

# Compile executable
printf "\n>>> Compiling MiniDATAMCNtuple_RunData2022B...\n"
COMPILER=$(root-config --cxx)
FLAGS=$(root-config --cflags --libs)
$COMPILER -g -O3 -o MiniDATAMCNtuple_RunData2022B  MiniDATAMCNtuple_RunData2022B.C MiniDATAMCNtuple.C $FLAGS

printf "\n>>> Running MiniDATAMCNtuple_RunData2022B\n"
./MiniDATAMCNtuple_RunData2022B
