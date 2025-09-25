#!/bin/bash

# Compile executable
printf "\n>>> Compiling MiniNtuple_RunTTToS...\n"
COMPILER=$(root-config --cxx)
FLAGS=$(root-config --cflags --libs)
$COMPILER -g -O3 -o MiniNtuple_RunTTToS  MiniNtuple_RunTTToS.C MiniNtuple.C $FLAGS

printf "\n>>> Running MiniNtuple_RunTTToS\n"
./MiniNtuple_RunTTToS
