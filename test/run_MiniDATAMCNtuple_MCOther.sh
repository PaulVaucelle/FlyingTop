#!/bin/bash

# Compile executable
printf "\n>>> Compiling MiniDATAMCNtuple_RunOther...\n"
COMPILER=$(root-config --cxx)
FLAGS=$(root-config --cflags --libs)
$COMPILER -g -O3 -o MiniDATAMCNtuple_RunOther  MiniDATAMCNtuple_RunOther.C MiniDATAMCNtuple.C $FLAGS

printf "\n>>> Running MiniDATAMCNtuple_RunOther\n"
./MiniDATAMCNtuple_RunOther
