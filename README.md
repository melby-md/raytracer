# Ray Tracer

![Cornell Box rendering](cornell_box.bmp)

An offline path tracing-based renderer written in C++11.

## Building

If you are on Windows use [w64devkit](https://github.com/skeeto/w64devkit), on \*NIX install a c++ compiler.
Then run:

    ./build.sh

To change the compiler, set the `CXX` enviroment variable:

    CXX=clang++ ./build.sh

To create an optimized build:

    ./build.sh opt

## Running

Execute `raytracer` specifying the scene file:

    ./raytracer scenes/cornell_box.txt

The render will be in the `image.bmp` file. Other examples are available in the `scenes` folder.
