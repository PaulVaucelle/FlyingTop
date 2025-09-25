#!/bin/bash

# Compile executable
printf "\n>>> Compiling MiniNtuple_RunOther...\n"
COMPILER=$(root-config --cxx)
FLAGS=$(root-config --cflags --libs)
$COMPILER -g -O3 -o MiniNtuple_RunOther  MiniNtuple_RunOther.C MiniNtuple.C $FLAGS

printf "\n>>> Running MiniNtuple_RunOther\n"
./MiniNtuple_RunOther
