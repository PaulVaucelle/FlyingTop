#!/bin/bash

# Compile executable
printf "\n>>> Compiling MiniNtuple_RunDY...\n"
COMPILER=$(root-config --cxx)
FLAGS=$(root-config --cflags --libs)
$COMPILER -g -O3 -o MiniNtuple_RunDY  MiniNtuple_RunDY.C MiniNtuple.C $FLAGS

printf "\n>>> Running MiniNtuple_RunDY\n"
./MiniNtuple_RunDY
