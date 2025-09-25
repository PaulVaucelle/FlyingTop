#!/bin/bash

# Compile executable
printf "\n>>> Compiling MiniDATAMCNtuple_RunDY...\n"
COMPILER=$(root-config --cxx)
FLAGS=$(root-config --cflags --libs)
$COMPILER -g -O3 -o MiniDATAMCNtuple_RunDY  MiniDATAMCNtuple_RunDY.C MiniDATAMCNtuple.C $FLAGS

printf "\n>>> Running MiniDATAMCNtuple_RunDY\n"
./MiniDATAMCNtuple_RunDY
