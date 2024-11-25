# LDPC Decoder Project

## Overview
This project implements an LDPC (Low-Density Parity-Check) decoder using C++. It includes functionality for generating LDPC codes, decoding with different LP relaxations.

- **src/**: Source codes.
- **length_120/**: Parity check matrix of code length 120 used in paper.
- **length_60/**: Parity check matrix of code length 60 used in paper.

## Dependencies
- C++ compiler supporting C++14 or higher.
- [Gurobi Optimizer](https://www.gurobi.com).

## Build the project:
Change the directory to 'src':
```bash
cd src
```
In `CMakeLists.txt`, update `GUROBI_HOME`:
```cmake
set(GUROBI_HOME "/path/to/your/gurobi/installation")
```
Then, run the following command to build the project:
```bash
mkdir build && cd build
cmake ..
make
```

## Execution:
Make a directory to store results:
```bash
mkdir output 
```

Execute with required parameters:
```bash
./LDPC_Decoder <clique_size> <density> <code_length> <num_trials>
```

Example:
```bash
./LDPC_Decoder 4 4 30 400
```

Parameters:
- clique_size: Size of cliques in LDPC graph
- density: Density parameter for LDPC code
- code_length: Length of the code
- num_trials: Number of trials to run

Output files will be generated in the `output/` directory.