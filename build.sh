#!/bin/sh -e
: "${CXX:=c++}"
set -x
${CXX} -o raytracer main.cpp -pipe -std=c++11 -fno-exceptions -fno-rtti -fopenmp -g3 -pedantic -Wall -Wextra -Wshadow -Wno-unused-parameter -Wno-unused-function -Wconversion -Wno-sign-conversion -Wdouble-promotion -lm
