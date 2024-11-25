#pragma once
#include <vector>
#include <unordered_map>
#include <string>

class generator {
public:
	generator(int num_bit);
	~generator();

	std::unordered_map<std::string, int> get_mapping();
	std::vector<std::vector<int>> get_check_matrix();
	std::vector<int> get_code();
	std::vector<int> get_noisy_code();

	void generate_LDPC(int clique_size, int num_density, int num_block, std::string path);
	void generate_CODE(int num_bit, int num_trials, std::string path);
	void read_LDPC(std::string filename, int clique_size, int num_density, int num_block);
	void read_code(std::string filename, int num_bit);

private:
	std::vector<std::vector<int>> check_matrix;
	std::vector<int> original_code;
	std::vector<int> received_code;
	std::unordered_map<std::string, int> mapping;
};