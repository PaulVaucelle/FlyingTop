#!/bin/bash

# Compile executable
printf "\n>>> Compiling MiniNtuple_RunTTTo2L2Nu...\n"
COMPILER=$(root-config --cxx)
FLAGS=$(root-config --cflags --libs)
$COMPILER -g -O3 -o MiniNtuple_RunTTTo2L2Nu  MiniNtuple_RunTTTo2L2Nu.C MiniNtuple.C $FLAGS

printf "\n>>> Running MiniNtuple_RunTTTo2L2Nu\n"
./MiniNtuple_RunTTTo2L2Nu
