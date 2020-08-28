# GridKit™

This is experimental code for prototyping interfaces for dynamic 
simulations, sensitivity analysis and optimization. Target applications 
are power grids, but the methodology and the framework could be used 
in other areas without major modifications. 

## Installation Guide

This guide is intended to give a guideline for installing GridKit™ and it's dependencies.

### Dependencies
You should have all of the following installed before installing GridKit™
- A version of
	- Ipopt >= 3.x
	- Sundials >= 4.x
	- Suitesparse >= 5.x
- CMake >= 3.12
- C++ 11 compiler

### Installing
It is recommended that you build and install GridKit™ in seperate directories in order to keep the main working tree clean, and so we suggest using a file structure like the following:
```
project_dir
└───source_dir  
│   │   ComponenetLib
│   │   Examples
│   │   Solver
|        ...
└───build_dir
└───install_dir
```
While there are many ways you can go about building the project, using the current example you could use the following set of commands:
```bash
cd project_dir
cmake -S source_dir -B build_dir
ccmake ./build_dir # In order to bring up the CMake curses frontend
cmake-gui ./build_dir # In order to bring up the CMake GUI frontend
cmake --build build_dir # To build all targets
cd build_dir; ctest # To validate examples were built successfully
cmake --build build_dir --target install_dir # To install - this is currently broken
```
Since the installation of Ipopt, Sundials and Suitesparse are all outside of the project_dir, you will need to specify the install directories of these modules within the CMake GUI. Once the root directory of these modules is specified in the CMake GUI, the include directories and the associated libraries should automatically be found. If you can load any modules beforehand, they will automatically be populated and you will not need to specify the root directory.

### Testing

Several examples are built when running `cmake --build buid_dir`, and these tests are automatically executed by running `ctest` in the build directory. The executables for each test can be found in `<build_dir>/Examples/<TestName>/<testexecutable>`.
