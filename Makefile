CXX = c++
WARNINGS = -pedantic -Wall -Wextra -Wconversion -Wdouble-promotion -Wshadow -Wno-unused-parameter -Wno-unused-function -Wno-sign-conversion -Werror
FLAGS = -pipe
OPT = -O3
CFLAGS = $(OPT) -std=c++11 -g3 -fno-exceptions -fno-rtti $(WARNINGS)
LDFLAGS = -lm

all: raytracer

raytracer: main.cpp common.hpp mathlib.hpp bmp.hpp
	$(CXX) -o $@ $< $(FLAGS) $(CFLAGS) $(LDFLAGS)

clean:
	rm -rf raytracer

.PHONY: all clean
