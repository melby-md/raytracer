# Ray Tracer

An offline path tracing-based renderer written in C++11.

## Building

If you are on Windows use [w64devkit](https://github.com/skeeto/w64devkit), on \*NIX install a c++ compiler.
Then run:

    ./build.sh

To change the compiler, set the `CXX` enviroment variable:

    CXX=clang++ ./build.sh

## Running

Execute on the same folder where the file `scene.txt` is:

    ./raytracer

To change the render modify `scene.txt`.
