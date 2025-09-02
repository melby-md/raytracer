#!/bin/sh -e
: "${CXX:=c++}"
set -x
${CXX} -o raytracer main.cpp -Og -pipe -std=c++11 -fno-exceptions -fno-rtti -fopenmp -g3 -Wpedantic -Wall -Wextra -Wshadow -Wconversion -Wdeprecated -Wdouble-promotion -Wno-unused-parameter -Wno-unused-function -Wno-sign-conversion -lm
