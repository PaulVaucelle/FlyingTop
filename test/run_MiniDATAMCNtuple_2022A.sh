#!/bin/bash

# Compile executable
printf "\n>>> Compiling MiniDATAMCNtuple_2022...\n"
COMPILER=$(root-config --cxx)
FLAGS=$(root-config --cflags --libs)
$COMPILER -g -O3 -o MiniDATAMCNtuple_RunData2022A  MiniDATAMCNtuple_RunData2022A.C MiniDATAMCNtuple.C $FLAGS

printf "\n>>> Running MiniDATAMCNtuple_S2022A\n"
./MiniDATAMCNtuple_RunData2022A
