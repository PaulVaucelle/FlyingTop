#!/bin/bash

# Compile executable
printf "\n>>> Compiling MiniNtuple_2022...\n"
COMPILER=$(root-config --cxx)
FLAGS=$(root-config --cflags --libs)
$COMPILER -g -O3 -o MiniNtuple_2022  MiniNtuple_RunS.C MiniNtuple.C $FLAGS

printf "\n>>> Running MiniNtuple_2022\n"
./MiniNtuple_2022
