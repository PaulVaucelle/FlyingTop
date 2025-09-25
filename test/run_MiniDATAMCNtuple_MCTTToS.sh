#!/bin/bash

# Compile executable
printf "\n>>> Compiling MiniDATAMCNtuple_RunTTToS...\n"
COMPILER=$(root-config --cxx)
FLAGS=$(root-config --cflags --libs)
$COMPILER -g -O3 -o MiniDATAMCNtuple_RunTTToS  MiniDATAMCNtuple_RunTTToS.C MiniDATAMCNtuple.C $FLAGS

printf "\n>>> Running MiniDATAMCNtuple_RunTTToS\n"
./MiniDATAMCNtuple_RunTTToS
