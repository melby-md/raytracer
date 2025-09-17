#!/bin/sh -e
cd "$(dirname "$0")"
: "${CXX:=c++}"
OPT=
for arg; do
	case "$arg" in
	opt) OPT='-O3 -DOPT';;
	esac
done
set -x
$CXX -o raytracer main.cpp $OPT -pipe -std=c++11 -fno-exceptions -fno-rtti \
     -fopenmp -g3 -Wpedantic -Wall -Wextra -Wshadow -Wconversion -Wdeprecated \
     -Wdouble-promotion -Wno-unused-parameter -Wno-unused-function \
     -Wno-sign-conversion -Wno-missing-braces -lm
