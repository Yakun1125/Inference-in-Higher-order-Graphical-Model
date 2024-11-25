#include "Generator.h"
#include "Optimizer.h"
#include <iostream>

int main(int argc, char* argv[]) {
    if (argc != 5) {
        std::cout << "Usage: " << argv[0] << " <clique_size> <density> <code_length> <num_trials>\n";
        return 1;
    }

    int clique_size = std::stoi(argv[1]);
    int density = std::stoi(argv[2]);
    int code_length = std::stoi(argv[3]);
    int num_trials = std::stoi(argv[4]);

    generator* input_code = new generator(clique_size * code_length);
    input_code->generate_LDPC(clique_size, density, code_length, "output/");
    input_code->generate_CODE(clique_size * code_length, num_trials, "output/");
    
    std::string matrix_file = "output/check_matrix_" + std::to_string(clique_size) + "_" + 
                             std::to_string(density) + "_" + std::to_string(code_length) + ".txt";
    input_code->read_LDPC(matrix_file, clique_size, density-1, code_length);
    std::cout<< " Code generated "<< std::endl;

    
    decoder* dec = new decoder(input_code, "output/");
    dec->parLP(input_code, 0, "output/Parity_IP.txt", "LP", num_trials);

    delete dec;
    delete input_code;
    return 0;
}