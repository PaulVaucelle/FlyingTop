#!/bin/bash

# Compile executable
printf "\n>>> Compiling MiniSecIntNtuple_2023...\n"
COMPILER=$(root-config --cxx)
FLAGS=$(root-config --cflags --libs)
$COMPILER -g -O3 -o MiniSecIntNtuple_2023  MiniSecIntNtuple_RunData2023.C MiniSecIntNtuple.C $FLAGS

printf "\n>>> Running MiniSecIntNtuple_2023\n"
./MiniSecIntNtuple_2023
