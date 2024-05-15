#include "Generator.h"
#include <random>
#include <fstream>
#include <iostream>

generator::generator(int num_bit) {
	// intialize everything
	//original code is all zero, and has length clique_size*num_block
	original_code.resize(num_bit, 0);
	received_code.resize(num_bit);
}

generator::~generator() {
	// destructor
}

std::vector<int> generator::get_code() {
	return original_code;
}

std::vector<int> generator::get_noisy_code() {
	return received_code;
}

std::vector<std::vector<int>> generator::get_check_matrix() {
	return check_matrix;
}

std::unordered_map<std::string, int> generator::get_mapping() {
	return mapping;
}

void generator::generate_CODE(int num_bit, std::string path) {
	std::string filename;
	// generate codeword sequence
	for (int co = 0; co < 31; co++) {
		double corrupt = 0.00 + co * 0.01;
		for (int num_sol = 0; num_sol < 400; num_sol++) {
			filename = path;
			filename += "code_" + std::to_string(num_bit) + "_" + std::to_string(corrupt) + "_" + std::to_string(num_sol) + ".txt";


			std::random_device rd2;
			std::mt19937 gen2(rd2());
			std::uniform_real_distribution<double> distribution(0, 1);

			std::vector<int> received_code(num_bit, 0);

			for (int i = 0; i < num_bit; i++) {
				double randomNumber = distribution(gen2);
				if (randomNumber < corrupt) {
					received_code[i] = -received_code[i] + 1;
				}
			}

			std::ofstream f(filename);
			for (int i = 0; i < num_bit; i++) {
				f << received_code[i] << " ";
			}
			f.close();

		}
	}
}

void generator::generate_LDPC(int clique_size, int num_density, int num_block, std::string path) {
	int num_bit = clique_size*num_block;
	int num_clique = num_density*num_block;

	/* Randomly Generate Parity Check Matrix
	* First generate sub matrix that contain cliques_size 1's in each row
	* Then permute the columns, generate num_density-1 sub matrices.
	* Every column would has num_density 1.
	*/
	std::vector<std::vector<int>> subm1(num_block, std::vector<int>(num_bit, 0));

	for (int row = 0; row < num_block; ++row) {
		for (int col = (row * clique_size); col < (row + 1) * clique_size; ++col) {
			subm1[row][col] = 1;
		}
	}

	check_matrix = subm1;

	std::random_device rd;
	std::mt19937 gen(rd());

	for (int group = 1; group < num_density; ++group) {
		std::vector<int> r(num_bit);
		for (int i = 0; i < num_bit; ++i) {
			r[i] = i;
		}
		std::shuffle(r.begin(), r.end(), gen);

		std::vector<std::vector<int>> temp_matrix(num_block, std::vector<int>(num_bit, 0));
		for (int row = 0; row < num_block; ++row) {
			for (int col = 0; col < num_bit; ++col) {
				temp_matrix[row][col] = subm1[row][r[col]];
			}
		}
		check_matrix.insert(check_matrix.end(), temp_matrix.begin(), temp_matrix.end());
	}

	// Save the parity check matrix to a file
	// named with num_density, num_block, clique_size

	std::string filename = path;
	filename += "check_matrix_"+std::to_string(clique_size) + "_" + std::to_string(num_density) + "_" + std::to_string(num_block) + ".txt";
	std::ofstream f(filename);

	for (int i = 0; i < num_clique; ++i) {
		for (int j = 0; j < num_bit; ++j) {
			f << check_matrix[i][j] << " ";
		}
		f << std::endl;
	}
	f.close();
}

void generator::read_LDPC(std::string filename, int clique_size, int num_density, int num_block) {
	std::ifstream fin(filename);
	// if the file is not found, return
	if (!fin) {
		return;
	}

	int num_bit = clique_size*num_block;
	int num_clique = num_density*num_block;

	check_matrix = std::vector<std::vector<int>>(num_clique, std::vector<int>(num_bit, 0));
	mapping.clear();

	for (int i = 0; i < num_clique; ++i) {
		for (int j = 0; j < num_bit; ++j) {
			fin >> check_matrix[i][j];
		}
	}

	fin.close();

	int _edge_idx = 0;
	for (const auto& row : check_matrix) {
		std::vector<int> arr;
		for (int i = 0; i < row.size(); i++) {
			if (row[i] == 1) {
				arr.push_back(i);
			}
		}
		std::sort(arr.begin(), arr.end());
		for (int size = 2; size <= clique_size; ++size) {
			std::vector<bool> v(arr.size());
			std::fill(v.end() - size, v.end(), true);

			do {
				std::string key;
				for (int i = 0; i < arr.size(); ++i) {
					if (v[i]) {
						key += std::to_string(arr[i]) + ",";
					}
				}
				if (mapping.find(key) == mapping.end()) {
					mapping[key] = _edge_idx++;
				}
			} while (std::next_permutation(v.begin(), v.end()));
		}
	}
}

void generator::read_code(std::string filename, int num_bit) {
	std::ifstream fin(filename);
	if (!fin) {
		return;
	}

	received_code = std::vector<int>(num_bit, 0);

	for (int i = 0; i < num_bit; i++) {
		fin >> received_code[i];
	}
}