#!/bin/sh -e
: "${CXX:=c++}"
FLAGS='-Og'
if [ "$1" = 'release' ]; then
	FLAGS='-O3 -DRELEASE'
fi
set -x
${CXX} -o raytracer main.cpp $FLAGS -pipe -std=c++11 -fno-exceptions -fno-rtti -fopenmp -g3 -Wpedantic -Wall -Wextra -Wshadow -Wconversion -Wdeprecated -Wdouble-promotion -Wno-unused-parameter -Wno-unused-function -Wno-sign-conversion -lm
