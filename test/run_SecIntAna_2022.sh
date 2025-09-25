#!/bin/bash

# Compile executable
printf "\n>>> Compiling MiniSecIntNtuple_2022...\n"
COMPILER=$(root-config --cxx)
FLAGS=$(root-config --cflags --libs)
$COMPILER -g -O3 -o MiniSecIntNtuple_2022  MiniSecIntNtuple_RunData2022.C MiniSecIntNtuple.C $FLAGS

printf "\n>>> Running MiniSecIntNtuple_2022\n"
./MiniSecIntNtuple_2022
