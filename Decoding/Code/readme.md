# LDPC Decoder Project

## Overview
This project implements an LDPC (Low-Density Parity-Check) decoder using C++. It includes functionality for generating LDPC codes, decoding with different LP relaxations.

## Dependencies
- C++ compiler supporting C++11 or higher.
- [Gurobi Optimizer](https://www.gurobi.com) - A powerful solver for linear programming and mixed-integer programming.


## Example Usage
Here is how you can generate and decode LDPC codes:
```cpp
generator* input_code = new generator(120);
input_code->generate_LDPC(4, 4, 30, "output/");
input_code->generate_CODE(120, "output/");
input_code->read_LDPC(("output/check_matrix_4_4_30.txt", 4, 3, 30);)
decoder* dec = new decoder(input_code, "output/");
dec->parLP(input_code, 0, "output/43output_Parity_IP.txt", "LP", 400);
// Cleanup
delete dec;
delete input_code;