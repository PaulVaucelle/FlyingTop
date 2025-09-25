#!/bin/bash

# Compile executable
printf "\n>>> Compiling MiniDATAMCNtuple_RunTTTo2L2Nu...\n"
COMPILER=$(root-config --cxx)
FLAGS=$(root-config --cflags --libs)
$COMPILER -g -O3 -o MiniDATAMCNtuple_RunTTTo2L2Nu  MiniDATAMCNtuple_RunTTTo2L2Nu.C MiniDATAMCNtuple.C $FLAGS

printf "\n>>> Running MiniDATAMCNtuple_RunTTTo2L2Nu\n"
./MiniDATAMCNtuple_RunTTTo2L2Nu
