CXX = c++
WARNINGS = -pedantic -Wall -Wextra -Wconversion -Wdouble-promotion -Wshadow -Wno-unused-parameter -Wno-unused-function -Wno-sign-conversion -Werror
FLAGS = -pipe
CFLAGS = -std=c++11 -g3 -fwhole-program -O3 -fno-exceptions -fno-rtti $(WARNINGS)
LDFLAGS = -lm

all: raytracer

raytracer: main.cpp common.hpp vector.hpp bmp.hpp
	$(CXX) -o $@ $< $(FLAGS) $(CFLAGS) $(LDFLAGS)

clean:
	rm -rf raytracer

.PHONY: all clean
