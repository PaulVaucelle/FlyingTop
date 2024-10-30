#!/bin/bash

# Compile executable
printf "\n>>> Compiling MiniSecIntNtuple_2024...\n"
COMPILER=$(root-config --cxx)
FLAGS=$(root-config --cflags --libs)
$COMPILER -g -O3 -o MiniSecIntNtuple_2024  MiniSecIntNtuple_RunData2024.C MiniSecIntNtuple.C $FLAGS

printf "\n>>> Running MiniSecIntNtuple_2024\n"
./MiniSecIntNtuple_2024
