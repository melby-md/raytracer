#!/bin/sh -e
: "${CXX:=c++}"
if [ $1 = opt ]; then
	FLAGS='-O3 -DOPT'
fi
set -x
${CXX} -o raytracer main.cpp $FLAGS -pipe -std=c++11 -fno-exceptions -fno-rtti -fopenmp -g3 -Wpedantic -Wall -Wextra -Wshadow -Wconversion -Wdeprecated -Wdouble-promotion -Wno-unused-parameter -Wno-unused-function -Wno-sign-conversion -Wno-missing-braces -lm
