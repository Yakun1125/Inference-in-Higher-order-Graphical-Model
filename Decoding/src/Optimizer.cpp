#include "Optimizer.h"
#include <algorithm>
#include <unordered_set>
#include "gurobi_c++.h"
#include <fstream>
#include <random>
#include <iomanip>

decoder::decoder(generator* input_code, std::string filename) {
	mapping = input_code->get_mapping();
	this->filename = filename;
}

int decoder::e(int i, int j, int k, int l) {
	std::vector<int> arr;
	if (i != -1) { arr.push_back(i); }
	if (j != -1) { arr.push_back(j); }
	if (k != -1) { arr.push_back(k); }
	if (l != -1) { arr.push_back(l); }
	std::sort(arr.begin(), arr.end());

	std::string key;
	for (int value : arr) {
		key += std::to_string(value) + ",";
	}

	return mapping[key];
}

int decoder::edge(std::vector<int> arr) {
	std::sort(arr.begin(), arr.end());

	std::string key;
	for (int value : arr) {
		key += std::to_string(value) + ",";
	}

	return mapping[key];
}

bool decoder::findCycle(std::vector<std::vector<int>>& sets, std::vector<int>& visited_sets, std::vector<int>& cycle, int startSet, int currentSet) {
	if (cycle.size() == sets.size() - 1) {
		// Check if the last set shares an element with the starting set
		for (int elem : sets[currentSet]) {
			if (std::find(sets[startSet].begin(), sets[startSet].end(), elem) != sets[startSet].end() && std::find(cycle.begin(), cycle.end(), elem) == cycle.end()) {
				cycle.push_back(elem); // Complete the cycle
				return true;
			}
		}
		return false; // No valid cycle found
	}

	for (int element : sets[currentSet]) {
		if (element != -1) { // Element is not used yet
			int nextSet = -1;
			int nextElem = -1;
			for (int i = 0; i < sets.size(); ++i) {// iterate through other sets, find if element exists
				auto next = std::find(sets[i].begin(), sets[i].end(), element);
				if (i != currentSet && visited_sets[i] != 1 && next != sets[i].end()) {
					nextElem = next - sets[i].begin();
					nextSet = i;
					break;
				}
			}

			if (nextSet != -1) {
				cycle.push_back(element);
				int savedElement = element;
				int savedNextElem = sets[nextSet][nextElem];
				element = -1; // Mark as used
				sets[nextSet][nextElem] = -1;
				visited_sets[nextSet] = 1;// Mark as visited
				if (findCycle(sets, visited_sets, cycle, startSet, nextSet)) {
					return true;
				}
				visited_sets[nextSet] = -1; // Unmark
				element = savedElement; // Unmark
				sets[nextSet][nextElem] = savedNextElem;
				cycle.pop_back();
			}
		}
	}

	return false; // No valid cycle found
}

bool decoder::Check_Structure(std::vector<std::vector<int>> Shared_HyperGraphs, std::vector<int>& cycle) {

	std::vector<int> visited_sets(Shared_HyperGraphs.size(), -1);
	bool cycleFound = false;

	for (int i = 0; i < Shared_HyperGraphs.size() - 1; ++i) {
		cycle.clear();
		visited_sets[i] = 1;
		if (findCycle(Shared_HyperGraphs, visited_sets, cycle, i, i)) {
			cycleFound = true;
			break;
		}
	}

	return cycleFound;
}

bool decoder::findCommonElementAndExclude(const std::vector<std::vector<int>>& vectors, int& commonElement, std::vector<std::vector<int>>& modifiedVectors) {
	// Create a set for each element in the first vector
	std::unordered_set<int> commonSet(vectors[0].begin(), vectors[0].end());

	// Iterate through the other vectors and find the common element
	for (const std::vector<int>& vec : vectors) {
		std::unordered_set<int> currentSet(vec.begin(), vec.end());

		// Intersect the current set with the commonSet
		for (auto it = commonSet.begin(); it != commonSet.end(); ) {
			if (currentSet.count(*it) == 0) {
				it = commonSet.erase(it);
			}
			else {
				++it;
			}
		}
	}

	if (commonSet.empty()) {
		return false; // No common element found
	}

	commonElement = *commonSet.begin();
	modifiedVectors = vectors;

	// Exclude the common element from each vector
	for (std::vector<int>& vec : modifiedVectors) {
		vec.erase(remove(vec.begin(), vec.end(), commonElement), vec.end());
	}

	return true;
}

void decoder::generateCombinations(std::vector<int>& currentCombination, std::vector<std::vector<int>>& combinations, int n, int k, int start) {
	if (k == 0) {
		combinations.push_back(currentCombination);
		return;
	}

	for (int i = start; i < n; i++) {
		currentCombination.push_back(i);
		generateCombinations(currentCombination, combinations, n, k - 1, i + 1);
		currentCombination.pop_back();
	}
}

void decoder::stdLP(generator* input_code, int print_level, std::string output_file, std::string version, int num_trials) {
	int num_graph = input_code->get_check_matrix().size();
	int num_bit = input_code->get_check_matrix()[0].size();

	GRBEnv env = GRBEnv(true);
	GRBVar* zv = 0;
	GRBVar* ze = 0;
	env.set(GRB_IntParam_OutputFlag, 0);
	env.start();
	GRBModel model = GRBModel(env);
	model.set(GRB_IntParam_OutputFlag, 0);
	model.set(GRB_IntParam_LogToConsole, 0);

	zv = new GRBVar[num_bit];
	ze = new GRBVar[mapping.size()];


	if (version == "LP") {
		for (int i = 0; i < num_bit; i++) {
			zv[i] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
		}
		for (int i = 0; i < mapping.size(); i++) {
			ze[i] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
		}
	}
	else {
		for (int i = 0; i < num_bit; i++) {
			zv[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		}
		for (int i = 0; i < mapping.size(); i++) {
			ze[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		}
	}

	for (const auto& row : input_code->get_check_matrix()) {
		std::vector<int> arr;
		for (int i = 0; i < row.size(); i++) {
			if (row[i] == 1) {
				arr.push_back(i);
			}
		}
		int v1 = arr[0]; int v2 = arr[1]; int v3 = arr[2]; int v4 = arr[3];
		model.addConstr(ze[e(v1, v2)] <= zv[v1]); model.addConstr(ze[e(v1, v2)] <= zv[v2]); model.addConstr(ze[e(v1, v2)] >= zv[v1] + zv[v2] - 1);
		model.addConstr(ze[e(v1, v3)] <= zv[v1]); model.addConstr(ze[e(v1, v3)] <= zv[v3]); model.addConstr(ze[e(v1, v3)] >= zv[v1] + zv[v3] - 1);
		model.addConstr(ze[e(v1, v4)] <= zv[v1]); model.addConstr(ze[e(v1, v4)] <= zv[v4]); model.addConstr(ze[e(v1, v4)] >= zv[v1] + zv[v4] - 1);
		model.addConstr(ze[e(v2, v3)] <= zv[v2]); model.addConstr(ze[e(v2, v3)] <= zv[v3]); model.addConstr(ze[e(v2, v3)] >= zv[v2] + zv[v3] - 1);
		model.addConstr(ze[e(v2, v4)] <= zv[v2]); model.addConstr(ze[e(v2, v4)] <= zv[v4]); model.addConstr(ze[e(v2, v4)] >= zv[v2] + zv[v4] - 1);
		model.addConstr(ze[e(v3, v4)] <= zv[v3]); model.addConstr(ze[e(v3, v4)] <= zv[v4]); model.addConstr(ze[e(v3, v4)] >= zv[v3] + zv[v4] - 1);
		model.addConstr(ze[e(v1, v2, v3)] <= zv[v1]); model.addConstr(ze[e(v1, v2, v3)] <= zv[v2]); model.addConstr(ze[e(v1, v2, v3)] <= zv[v3]); model.addConstr(ze[e(v1, v2, v3)] >= zv[v1] + zv[v2] + zv[v3] - 2);
		model.addConstr(ze[e(v1, v2, v4)] <= zv[v1]); model.addConstr(ze[e(v1, v2, v4)] <= zv[v2]); model.addConstr(ze[e(v1, v2, v4)] <= zv[v4]); model.addConstr(ze[e(v1, v2, v4)] >= zv[v1] + zv[v2] + zv[v4] - 2);
		model.addConstr(ze[e(v1, v3, v4)] <= zv[v1]); model.addConstr(ze[e(v1, v3, v4)] <= zv[v3]); model.addConstr(ze[e(v1, v3, v4)] <= zv[v4]); model.addConstr(ze[e(v1, v3, v4)] >= zv[v1] + zv[v3] + zv[v4] - 2);
		model.addConstr(ze[e(v2, v3, v4)] <= zv[v2]); model.addConstr(ze[e(v2, v3, v4)] <= zv[v3]); model.addConstr(ze[e(v2, v3, v4)] <= zv[v4]); model.addConstr(ze[e(v2, v3, v4)] >= zv[v2] + zv[v3] + zv[v4] - 2);
		model.addConstr(ze[e(v1, v2, v3, v4)] <= zv[v1]); model.addConstr(ze[e(v1, v2, v3, v4)] <= zv[v2]); model.addConstr(ze[e(v1, v2, v3, v4)] <= zv[v3]); model.addConstr(ze[e(v1, v2, v3, v4)] <= zv[v4]); model.addConstr(ze[e(v1, v2, v3, v4)] >= zv[v1] + zv[v2] + zv[v3] + zv[v4] - 3);

		model.addConstr(-zv[v1] - zv[v2] - zv[v3] - zv[v4] + 2 * ze[e(v1, v2)] + 2 * ze[e(v1, v3)] + 2 * ze[e(v1, v4)] + 2 * ze[e(v2, v3)] + 2 * ze[e(v2, v4)] + 2 * ze[e(v3, v4)] - 4 * ze[e(v1, v2, v3)] - 4 * ze[e(v1, v2, v4)] - 4 * ze[e(v1, v3, v4)] - 4 * ze[e(v2, v3, v4)] + 8 * ze[e(v1, v2, v3, v4)] + 1 == 1);
	}

	std::ofstream fin(output_file);

	std::cout << "stdLP" << std::endl;

	for (int co = 0; co < 21; co++) {
		double corrupt = 0.00 + co * 0.01;
		std::cout << corrupt << std::endl;
		std::string temp = "corruption level: " + std::to_string(corrupt) + "\n";
		fin << temp;

		/* Record the tightness rate and recovery rate*/
		double tightness = 0;
		double partial_recovery = 0;
		double exact_recovery = 0;

		for (int num_sol = 0; num_sol < num_trials; num_sol++) {
			std::string path = filename;
			path += "code_" + std::to_string(num_bit) + "_" + std::to_string(corrupt) + "_" + std::to_string(num_sol) + ".txt";
			input_code->read_code(path, num_bit);
			GRBLinExpr obj = 1;

			for (int i = 0; i < num_bit; i++) {
				obj += (-2 * input_code->get_noisy_code()[i] + 1) * zv[i];
			}
			model.setObjective(obj, GRB_MINIMIZE);
			model.reset();
			auto begin = std::clock();
			model.optimize();
			auto end = std::clock();
			auto duration = double(end - begin) / CLOCKS_PER_SEC;
			// store results separated by space
			fin << " instance " << num_sol << ":";
			fin << duration << " ";
			double opt = obj.getValue();
			fin << opt << " ";

			double recovery_rate = 0;
			bool is_recovered = true;
			bool is_binary = true;
			for (int i = 0; i < num_bit; i++) {
				if (zv[i].get(GRB_DoubleAttr_X) > 0.0001 && zv[i].get(GRB_DoubleAttr_X) < 0.9999) {
					is_binary = false;
					is_recovered = false;
					break;
				}
			}
			if (is_binary) {
				for (int i = 0; i < num_bit; i++) {
					if (zv[i].get(GRB_DoubleAttr_X) < 0.0001) {
						recovery_rate += 1;
					}
					else {
						is_recovered = false;
					}
				}
			}
			else {
				for (int i = 0; i < num_bit; i++) {
					if (std::round(zv[i].get(GRB_DoubleAttr_X)) < 0.0001) {
						recovery_rate += 1;
					}
				}
			}

			fin << is_binary << " ";
			fin << recovery_rate / num_bit << " ";
			fin << is_recovered << " ";
			fin << "\n";

			tightness += (int)is_binary;
			exact_recovery += is_recovered;
			partial_recovery += recovery_rate;

			if (print_level >= 1) {
				std::cout << " clique: ";
				std::cout << " time: " << duration << " obj value: " << obj.getValue() << " gap: " << std::abs(obj.getValue() - opt) / std::abs(opt) << std::endl;
				if (print_level >= 2) {
					std::cout << " Coding Result: ";
					for (int i = 0; i < num_bit; i++) {
						std::cout << zv[i].get(GRB_DoubleAttr_X) << " ";
					}
					std::cout << std::endl;
				}
				std::cout << "   Is solution recovered: " << is_recovered << std::endl;

				if (print_level >= 3) {

					std::cout << " Received Code: ";
					for (int i = 0; i < num_bit; i++) {
						std::cout << input_code->get_noisy_code()[i] << " ";
					}
					std::cout << std::endl;

					std::cout << " Original Code: ";
					for (int i = 0; i < num_bit; i++) {
						std::cout << input_code->get_code()[i] << " ";
					}
					std::cout << std::endl;

					bool parity_check = true;
					for (const auto& row : input_code->get_check_matrix()) {
						int temp = 0;
						for (int i = 0; i < num_bit; i++) {
							temp += input_code->get_noisy_code()[i] * row[i];
						}
						if (temp % 2 != 0) {
							parity_check = false;
						}
					}

					std::cout << " Is noisy code satisfies the parity check: " << parity_check << std::endl;

					if (print_level >= 3) {
						std::cout << " Coding Result: ";
						for (int i = 0; i < num_bit; i++) {
							std::cout << zv[i].get(GRB_DoubleAttr_X) << " ";
						}
						std::cout << std::endl;
					}
				}
			}
		}

		std::cout << "tightness rate " << tightness / num_trials << " ";
		std::cout << "partial recovery rate " << partial_recovery / (num_trials * num_bit) <<" ";
		std::cout << "exact recovery rate " << exact_recovery / num_trials << std::endl;
	}

	// clean up
	delete[] zv;
	delete[] ze;
}

void decoder::flLP(generator* input_code, int print_level, std::string output_file, std::string version, int num_trials) {
	int num_graph = input_code->get_check_matrix().size();
	int num_bit = input_code->get_check_matrix()[0].size();

	GRBEnv env = GRBEnv(true);
	GRBVar* zv = 0;
	GRBVar* ze = 0;
	env.set(GRB_IntParam_OutputFlag, 0);
	env.start();
	GRBModel model = GRBModel(env);
	model.set(GRB_IntParam_OutputFlag, 0);
	model.set(GRB_IntParam_LogToConsole, 0);

	zv = new GRBVar[num_bit];
	ze = new GRBVar[mapping.size()];


	if (version == "LP") {
		for (int i = 0; i < num_bit; i++) {
			zv[i] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
		}
		for (int i = 0; i < mapping.size(); i++) {
			ze[i] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
		}
	}
	else {
		for (int i = 0; i < num_bit; i++) {
			zv[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		}
		for (int i = 0; i < mapping.size(); i++) {
			ze[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		}
	}

	for (const auto& row : input_code->get_check_matrix()) {
		std::vector<int> arr;
		for (int i = 0; i < row.size(); i++) {
			if (row[i] == 1) {
				arr.push_back(i);
			}
		}
		int v1 = arr[0]; int v2 = arr[1]; int v3 = arr[2]; int v4 = arr[3];
		model.addConstr(ze[e(v1, v2)] <= zv[v1]); model.addConstr(ze[e(v1, v2)] <= zv[v2]); model.addConstr(ze[e(v1, v2)] >= zv[v1] + zv[v2] - 1);
		model.addConstr(ze[e(v1, v3)] <= zv[v1]); model.addConstr(ze[e(v1, v3)] <= zv[v3]); model.addConstr(ze[e(v1, v3)] >= zv[v1] + zv[v3] - 1);
		model.addConstr(ze[e(v1, v4)] <= zv[v1]); model.addConstr(ze[e(v1, v4)] <= zv[v4]); model.addConstr(ze[e(v1, v4)] >= zv[v1] + zv[v4] - 1);
		model.addConstr(ze[e(v2, v3)] <= zv[v2]); model.addConstr(ze[e(v2, v3)] <= zv[v3]); model.addConstr(ze[e(v2, v3)] >= zv[v2] + zv[v3] - 1);
		model.addConstr(ze[e(v2, v4)] <= zv[v2]); model.addConstr(ze[e(v2, v4)] <= zv[v4]); model.addConstr(ze[e(v2, v4)] >= zv[v2] + zv[v4] - 1);
		model.addConstr(ze[e(v3, v4)] <= zv[v3]); model.addConstr(ze[e(v3, v4)] <= zv[v4]); model.addConstr(ze[e(v3, v4)] >= zv[v3] + zv[v4] - 1);
		model.addConstr(ze[e(v1, v2, v3)] <= zv[v1]); model.addConstr(ze[e(v1, v2, v3)] <= zv[v2]); model.addConstr(ze[e(v1, v2, v3)] <= zv[v3]); model.addConstr(ze[e(v1, v2, v3)] >= zv[v1] + zv[v2] + zv[v3] - 2);
		model.addConstr(ze[e(v1, v2, v4)] <= zv[v1]); model.addConstr(ze[e(v1, v2, v4)] <= zv[v2]); model.addConstr(ze[e(v1, v2, v4)] <= zv[v4]); model.addConstr(ze[e(v1, v2, v4)] >= zv[v1] + zv[v2] + zv[v4] - 2);
		model.addConstr(ze[e(v1, v3, v4)] <= zv[v1]); model.addConstr(ze[e(v1, v3, v4)] <= zv[v3]); model.addConstr(ze[e(v1, v3, v4)] <= zv[v4]); model.addConstr(ze[e(v1, v3, v4)] >= zv[v1] + zv[v3] + zv[v4] - 2);
		model.addConstr(ze[e(v2, v3, v4)] <= zv[v2]); model.addConstr(ze[e(v2, v3, v4)] <= zv[v3]); model.addConstr(ze[e(v2, v3, v4)] <= zv[v4]); model.addConstr(ze[e(v2, v3, v4)] >= zv[v2] + zv[v3] + zv[v4] - 2);
		model.addConstr(ze[e(v1, v2, v3, v4)] <= zv[v1]); model.addConstr(ze[e(v1, v2, v3, v4)] <= zv[v2]); model.addConstr(ze[e(v1, v2, v3, v4)] <= zv[v3]); model.addConstr(ze[e(v1, v2, v3, v4)] <= zv[v4]); model.addConstr(ze[e(v1, v2, v3, v4)] >= zv[v1] + zv[v2] + zv[v3] + zv[v4] - 3);

		model.addConstr(ze[e(v1, v2)] <= ze[e(v1, v2, v3)] - zv[v3] + 1); model.addConstr(ze[e(v1, v2, v3)] <= ze[e(v1, v2)]);
		model.addConstr(ze[e(v1, v2)] <= ze[e(v1, v2, v4)] - zv[v4] + 1); model.addConstr(ze[e(v1, v2, v4)] <= ze[e(v1, v2)]);
		model.addConstr(ze[e(v1, v3)] <= ze[e(v1, v2, v3)] - zv[v2] + 1); model.addConstr(ze[e(v1, v2, v3)] <= ze[e(v1, v3)]);
		model.addConstr(ze[e(v1, v3)] <= ze[e(v1, v3, v4)] - zv[v4] + 1); model.addConstr(ze[e(v1, v3, v4)] <= ze[e(v1, v3)]);
		model.addConstr(ze[e(v1, v4)] <= ze[e(v1, v2, v4)] - zv[v2] + 1); model.addConstr(ze[e(v1, v2, v4)] <= ze[e(v1, v4)]);
		model.addConstr(ze[e(v1, v4)] <= ze[e(v1, v3, v4)] - zv[v3] + 1); model.addConstr(ze[e(v1, v3, v4)] <= ze[e(v1, v4)]);
		model.addConstr(ze[e(v2, v3)] <= ze[e(v1, v2, v3)] - zv[v1] + 1); model.addConstr(ze[e(v1, v2, v3)] <= ze[e(v2, v3)]);
		model.addConstr(ze[e(v2, v3)] <= ze[e(v2, v3, v4)] - zv[v4] + 1); model.addConstr(ze[e(v2, v3, v4)] <= ze[e(v2, v3)]);
		model.addConstr(ze[e(v2, v4)] <= ze[e(v1, v2, v4)] - zv[v1] + 1); model.addConstr(ze[e(v1, v2, v4)] <= ze[e(v2, v4)]);
		model.addConstr(ze[e(v2, v4)] <= ze[e(v2, v3, v4)] - zv[v3] + 1); model.addConstr(ze[e(v2, v3, v4)] <= ze[e(v2, v4)]);
		model.addConstr(ze[e(v3, v4)] <= ze[e(v1, v3, v4)] - zv[v1] + 1); model.addConstr(ze[e(v1, v3, v4)] <= ze[e(v3, v4)]);
		model.addConstr(ze[e(v3, v4)] <= ze[e(v2, v3, v4)] - zv[v2] + 1); model.addConstr(ze[e(v2, v3, v4)] <= ze[e(v3, v4)]);

		model.addConstr(ze[e(v1, v2, v3)] <= ze[e(v1, v2, v3, v4)] - zv[v4] + 1); model.addConstr(ze[e(v1, v2, v3, v4)] <= ze[e(v1, v2, v3)]);
		model.addConstr(ze[e(v1, v2, v4)] <= ze[e(v1, v2, v3, v4)] - zv[v3] + 1); model.addConstr(ze[e(v1, v2, v3, v4)] <= ze[e(v1, v2, v4)]);
		model.addConstr(ze[e(v1, v3, v4)] <= ze[e(v1, v2, v3, v4)] - zv[v2] + 1); model.addConstr(ze[e(v1, v2, v3, v4)] <= ze[e(v1, v3, v4)]);
		model.addConstr(ze[e(v2, v3, v4)] <= ze[e(v1, v2, v3, v4)] - zv[v1] + 1); model.addConstr(ze[e(v1, v2, v3, v4)] <= ze[e(v2, v3, v4)]);

		// flower
		model.addConstr(ze[e(v1, v2)] + ze[e(v3, v4)] - ze[e(v1, v2, v3, v4)] <= 1);
		model.addConstr(ze[e(v1, v3)] + ze[e(v2, v4)] - ze[e(v1, v2, v3, v4)] <= 1);
		model.addConstr(ze[e(v1, v4)] + ze[e(v2, v3)] - ze[e(v1, v2, v3, v4)] <= 1);

		model.addConstr(-zv[v1] - zv[v2] - zv[v3] - zv[v4] + 2 * ze[e(v1, v2)] + 2 * ze[e(v1, v3)] + 2 * ze[e(v1, v4)] + 2 * ze[e(v2, v3)] + 2 * ze[e(v2, v4)] + 2 * ze[e(v3, v4)] - 4 * ze[e(v1, v2, v3)] - 4 * ze[e(v1, v2, v4)] - 4 * ze[e(v1, v3, v4)] - 4 * ze[e(v2, v3, v4)] + 8 * ze[e(v1, v2, v3, v4)] + 1 == 1);
	}

	std::ofstream fin(output_file);

	std::cout << "flLP" << std::endl;

	for (int co = 0; co < 31; co++) {
		double corrupt = 0.00 + co * 0.01;
		std::cout << corrupt << std::endl;
		std::string temp = "corruption level: " + std::to_string(corrupt) + "\n";
		fin << temp;

		/* Record the tightness rate and recovery rate*/
		double tightness = 0;
		double partial_recovery = 0;
		double exact_recovery = 0;

		for (int num_sol = 0; num_sol < num_trials; num_sol++) {
			std::string path = filename;
			path += "code_" + std::to_string(num_bit) + "_" + std::to_string(corrupt) + "_" + std::to_string(num_sol) + ".txt";
			input_code->read_code(path, num_bit);
			GRBLinExpr obj = 1;

			for (int i = 0; i < num_bit; i++) {
				obj += (-2 * input_code->get_noisy_code()[i] + 1) * zv[i];
			}
			model.setObjective(obj, GRB_MINIMIZE);
			model.reset();
			auto begin = std::clock();
			model.optimize();
			auto end = std::clock();
			auto duration = double(end - begin) / CLOCKS_PER_SEC;
			// store results separated by space
			fin << " instance " << num_sol << ":";
			fin << duration << " ";
			double opt = obj.getValue();
			fin << opt << " ";

			double recovery_rate = 0;
			bool is_recovered = true;
			bool is_binary = true;
			for (int i = 0; i < num_bit; i++) {
				if (zv[i].get(GRB_DoubleAttr_X) > 0.0001 && zv[i].get(GRB_DoubleAttr_X) < 0.9999) {
					is_binary = false;
					is_recovered = false;
					break;
				}
			}
			if (is_binary) {
				for (int i = 0; i < num_bit; i++) {
					if (zv[i].get(GRB_DoubleAttr_X) < 0.0001) {
						recovery_rate += 1;
					}
					else {
						is_recovered = false;
					}
				}
			}
			else {
				for (int i = 0; i < num_bit; i++) {
					if (std::round(zv[i].get(GRB_DoubleAttr_X)) < 0.0001) {
						recovery_rate += 1;
					}
				}
			}

			fin << is_binary << " ";
			fin << recovery_rate / num_bit << " ";
			fin << is_recovered << " ";
			fin << "\n";

			tightness += (int)is_binary;
			exact_recovery += is_recovered;
			partial_recovery += recovery_rate;

			if (print_level >= 1) {
				std::cout << " clique: ";
				std::cout << " time: " << duration << " obj value: " << obj.getValue() << " gap: " << std::abs(obj.getValue() - opt) / std::abs(opt) << std::endl;
				if (print_level >= 2) {
					std::cout << " Coding Result: ";
					for (int i = 0; i < num_bit; i++) {
						std::cout << zv[i].get(GRB_DoubleAttr_X) << " ";
					}
					std::cout << std::endl;
				}
				std::cout << "   Is solution recovered: " << is_recovered << std::endl;

				if (print_level >= 3) {

					std::cout << " Received Code: ";
					for (int i = 0; i < num_bit; i++) {
						std::cout << input_code->get_noisy_code()[i] << " ";
					}
					std::cout << std::endl;

					std::cout << " Original Code: ";
					for (int i = 0; i < num_bit; i++) {
						std::cout << input_code->get_code()[i] << " ";
					}
					std::cout << std::endl;

					bool parity_check = true;
					for (const auto& row : input_code->get_check_matrix()) {
						int temp = 0;
						for (int i = 0; i < num_bit; i++) {
							temp += input_code->get_noisy_code()[i] * row[i];
						}
						if (temp % 2 != 0) {
							parity_check = false;
						}
					}

					std::cout << " Is noisy code satisfies the parity check: " << parity_check << std::endl;

					if (print_level >= 3) {
						std::cout << " Coding Result: ";
						for (int i = 0; i < num_bit; i++) {
							std::cout << zv[i].get(GRB_DoubleAttr_X) << " ";
						}
						std::cout << std::endl;
					}
				}
			}
		}

		std::cout << "tightness rate " << tightness / num_trials << " ";
		std::cout << "partial recovery rate " << partial_recovery / (num_trials * num_bit) << " ";
		std::cout << "exact recovery rate " << exact_recovery / num_trials << std::endl;
	}

	// clean up
	delete[] zv;
	delete[] ze;
}

void decoder::runLP(generator* input_code, int print_level, std::string output_file, std::string version, int num_trials) {
	int num_graph = input_code->get_check_matrix().size();
	int num_bit = input_code->get_check_matrix()[0].size();

	GRBEnv env = GRBEnv(true);
	GRBVar* zv = 0;
	GRBVar* ze = 0;
	env.set(GRB_IntParam_OutputFlag, 0);
	env.start();
	GRBModel model = GRBModel(env);
	model.set(GRB_IntParam_OutputFlag, 0);
	model.set(GRB_IntParam_LogToConsole, 0);

	zv = new GRBVar[num_bit];
	ze = new GRBVar[mapping.size()];


	if (version == "LP") {
		for (int i = 0; i < num_bit; i++) {
			zv[i] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
		}
		for (int i = 0; i < mapping.size(); i++) {
			ze[i] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
		}
	}
	else {
		for (int i = 0; i < num_bit; i++) {
			zv[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		}
		for (int i = 0; i < mapping.size(); i++) {
			ze[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		}
	}

	for (const auto& row : input_code->get_check_matrix()) {
		std::vector<int> arr;
		for (int i = 0; i < row.size(); i++) {
			if (row[i] == 1) {
				arr.push_back(i);
			}
		}
		int v1 = arr[0]; int v2 = arr[1]; int v3 = arr[2]; int v4 = arr[3];
		model.addConstr(ze[e(v1, v2)] <= zv[v1]); model.addConstr(ze[e(v1, v2)] <= zv[v2]); model.addConstr(ze[e(v1, v2)] >= zv[v1] + zv[v2] - 1);
		model.addConstr(ze[e(v1, v3)] <= zv[v1]); model.addConstr(ze[e(v1, v3)] <= zv[v3]); model.addConstr(ze[e(v1, v3)] >= zv[v1] + zv[v3] - 1);
		model.addConstr(ze[e(v1, v4)] <= zv[v1]); model.addConstr(ze[e(v1, v4)] <= zv[v4]); model.addConstr(ze[e(v1, v4)] >= zv[v1] + zv[v4] - 1);
		model.addConstr(ze[e(v2, v3)] <= zv[v2]); model.addConstr(ze[e(v2, v3)] <= zv[v3]); model.addConstr(ze[e(v2, v3)] >= zv[v2] + zv[v3] - 1);
		model.addConstr(ze[e(v2, v4)] <= zv[v2]); model.addConstr(ze[e(v2, v4)] <= zv[v4]); model.addConstr(ze[e(v2, v4)] >= zv[v2] + zv[v4] - 1);
		model.addConstr(ze[e(v3, v4)] <= zv[v3]); model.addConstr(ze[e(v3, v4)] <= zv[v4]); model.addConstr(ze[e(v3, v4)] >= zv[v3] + zv[v4] - 1);
		model.addConstr(ze[e(v1, v2, v3)] <= zv[v1]); model.addConstr(ze[e(v1, v2, v3)] <= zv[v2]); model.addConstr(ze[e(v1, v2, v3)] <= zv[v3]); model.addConstr(ze[e(v1, v2, v3)] >= zv[v1] + zv[v2] + zv[v3] - 2);
		model.addConstr(ze[e(v1, v2, v4)] <= zv[v1]); model.addConstr(ze[e(v1, v2, v4)] <= zv[v2]); model.addConstr(ze[e(v1, v2, v4)] <= zv[v4]); model.addConstr(ze[e(v1, v2, v4)] >= zv[v1] + zv[v2] + zv[v4] - 2);
		model.addConstr(ze[e(v1, v3, v4)] <= zv[v1]); model.addConstr(ze[e(v1, v3, v4)] <= zv[v3]); model.addConstr(ze[e(v1, v3, v4)] <= zv[v4]); model.addConstr(ze[e(v1, v3, v4)] >= zv[v1] + zv[v3] + zv[v4] - 2);
		model.addConstr(ze[e(v2, v3, v4)] <= zv[v2]); model.addConstr(ze[e(v2, v3, v4)] <= zv[v3]); model.addConstr(ze[e(v2, v3, v4)] <= zv[v4]); model.addConstr(ze[e(v2, v3, v4)] >= zv[v2] + zv[v3] + zv[v4] - 2);
		model.addConstr(ze[e(v1, v2, v3, v4)] <= zv[v1]); model.addConstr(ze[e(v1, v2, v3, v4)] <= zv[v2]); model.addConstr(ze[e(v1, v2, v3, v4)] <= zv[v3]); model.addConstr(ze[e(v1, v2, v3, v4)] <= zv[v4]); model.addConstr(ze[e(v1, v2, v3, v4)] >= zv[v1] + zv[v2] + zv[v3] + zv[v4] - 3);

		model.addConstr(ze[e(v1, v2)] <= ze[e(v1, v2, v3)] - zv[v3] + 1); model.addConstr(ze[e(v1, v2, v3)] <= ze[e(v1, v2)]);
		model.addConstr(ze[e(v1, v2)] <= ze[e(v1, v2, v4)] - zv[v4] + 1); model.addConstr(ze[e(v1, v2, v4)] <= ze[e(v1, v2)]);
		model.addConstr(ze[e(v1, v3)] <= ze[e(v1, v2, v3)] - zv[v2] + 1); model.addConstr(ze[e(v1, v2, v3)] <= ze[e(v1, v3)]);
		model.addConstr(ze[e(v1, v3)] <= ze[e(v1, v3, v4)] - zv[v4] + 1); model.addConstr(ze[e(v1, v3, v4)] <= ze[e(v1, v3)]);
		model.addConstr(ze[e(v1, v4)] <= ze[e(v1, v2, v4)] - zv[v2] + 1); model.addConstr(ze[e(v1, v2, v4)] <= ze[e(v1, v4)]);
		model.addConstr(ze[e(v1, v4)] <= ze[e(v1, v3, v4)] - zv[v3] + 1); model.addConstr(ze[e(v1, v3, v4)] <= ze[e(v1, v4)]);
		model.addConstr(ze[e(v2, v3)] <= ze[e(v1, v2, v3)] - zv[v1] + 1); model.addConstr(ze[e(v1, v2, v3)] <= ze[e(v2, v3)]);
		model.addConstr(ze[e(v2, v3)] <= ze[e(v2, v3, v4)] - zv[v4] + 1); model.addConstr(ze[e(v2, v3, v4)] <= ze[e(v2, v3)]);
		model.addConstr(ze[e(v2, v4)] <= ze[e(v1, v2, v4)] - zv[v1] + 1); model.addConstr(ze[e(v1, v2, v4)] <= ze[e(v2, v4)]);
		model.addConstr(ze[e(v2, v4)] <= ze[e(v2, v3, v4)] - zv[v3] + 1); model.addConstr(ze[e(v2, v3, v4)] <= ze[e(v2, v4)]);
		model.addConstr(ze[e(v3, v4)] <= ze[e(v1, v3, v4)] - zv[v1] + 1); model.addConstr(ze[e(v1, v3, v4)] <= ze[e(v3, v4)]);
		model.addConstr(ze[e(v3, v4)] <= ze[e(v2, v3, v4)] - zv[v2] + 1); model.addConstr(ze[e(v2, v3, v4)] <= ze[e(v3, v4)]);

		model.addConstr(ze[e(v1, v2, v3)] <= ze[e(v1, v2, v3, v4)] - zv[v4] + 1); model.addConstr(ze[e(v1, v2, v3, v4)] <= ze[e(v1, v2, v3)]);
		model.addConstr(ze[e(v1, v2, v4)] <= ze[e(v1, v2, v3, v4)] - zv[v3] + 1); model.addConstr(ze[e(v1, v2, v3, v4)] <= ze[e(v1, v2, v4)]);
		model.addConstr(ze[e(v1, v3, v4)] <= ze[e(v1, v2, v3, v4)] - zv[v2] + 1); model.addConstr(ze[e(v1, v2, v3, v4)] <= ze[e(v1, v3, v4)]);
		model.addConstr(ze[e(v2, v3, v4)] <= ze[e(v1, v2, v3, v4)] - zv[v1] + 1); model.addConstr(ze[e(v1, v2, v3, v4)] <= ze[e(v2, v3, v4)]);

		// flower
		model.addConstr(ze[e(v1, v2)] + ze[e(v3, v4)] - ze[e(v1, v2, v3, v4)] <= 1);
		model.addConstr(ze[e(v1, v3)] + ze[e(v2, v4)] - ze[e(v1, v2, v3, v4)] <= 1);
		model.addConstr(ze[e(v1, v4)] + ze[e(v2, v3)] - ze[e(v1, v2, v3, v4)] <= 1);

		// running intersection inequality
		model.addConstr(-zv[v1] + ze[e(v1, v2)] + ze[e(v1, v3)] - ze[e(v1, v2, v3)] <= 0);
		model.addConstr(-zv[v2] + ze[e(v1, v2)] + ze[e(v2, v3)] - ze[e(v1, v2, v3)] <= 0);
		model.addConstr(-zv[v3] + ze[e(v1, v3)] + ze[e(v2, v3)] - ze[e(v1, v2, v3)] <= 0);
		model.addConstr(-zv[v1] + ze[e(v1, v2)] + ze[e(v1, v4)] - ze[e(v1, v2, v4)] <= 0);
		model.addConstr(-zv[v2] + ze[e(v1, v2)] + ze[e(v2, v4)] - ze[e(v1, v2, v4)] <= 0);
		model.addConstr(-zv[v4] + ze[e(v1, v4)] + ze[e(v2, v4)] - ze[e(v1, v2, v4)] <= 0);
		model.addConstr(-zv[v1] + ze[e(v1, v3)] + ze[e(v1, v4)] - ze[e(v1, v3, v4)] <= 0);
		model.addConstr(-zv[v3] + ze[e(v1, v3)] + ze[e(v3, v4)] - ze[e(v1, v3, v4)] <= 0);
		model.addConstr(-zv[v4] + ze[e(v1, v4)] + ze[e(v3, v4)] - ze[e(v1, v3, v4)] <= 0);
		model.addConstr(-zv[v2] + ze[e(v2, v3)] + ze[e(v2, v4)] - ze[e(v2, v3, v4)] <= 0);
		model.addConstr(-zv[v3] + ze[e(v2, v3)] + ze[e(v3, v4)] - ze[e(v2, v3, v4)] <= 0);
		model.addConstr(-zv[v4] + ze[e(v2, v4)] + ze[e(v3, v4)] - ze[e(v2, v3, v4)] <= 0);

		model.addConstr(-zv[v1] + ze[e(v1, v2, v3)] + ze[e(v1, v2, v4)] - ze[e(v1, v2, v3, v4)] <= 0);
		model.addConstr(-zv[v2] + ze[e(v1, v2, v3)] + ze[e(v1, v2, v4)] - ze[e(v1, v2, v3, v4)] <= 0);
		model.addConstr(-zv[v1] + ze[e(v1, v2, v3)] + ze[e(v1, v3, v4)] - ze[e(v1, v2, v3, v4)] <= 0);
		model.addConstr(-zv[v3] + ze[e(v1, v2, v3)] + ze[e(v1, v3, v4)] - ze[e(v1, v2, v3, v4)] <= 0);
		model.addConstr(-zv[v2] + ze[e(v1, v2, v3)] + ze[e(v2, v3, v4)] - ze[e(v1, v2, v3, v4)] <= 0);
		model.addConstr(-zv[v3] + ze[e(v1, v2, v3)] + ze[e(v2, v3, v4)] - ze[e(v1, v2, v3, v4)] <= 0);
		model.addConstr(-zv[v1] + ze[e(v1, v2, v4)] + ze[e(v1, v3, v4)] - ze[e(v1, v2, v3, v4)] <= 0);
		model.addConstr(-zv[v4] + ze[e(v1, v2, v4)] + ze[e(v1, v3, v4)] - ze[e(v1, v2, v3, v4)] <= 0);
		model.addConstr(-zv[v2] + ze[e(v1, v2, v4)] + ze[e(v2, v3, v4)] - ze[e(v1, v2, v3, v4)] <= 0);
		model.addConstr(-zv[v4] + ze[e(v1, v2, v4)] + ze[e(v2, v3, v4)] - ze[e(v1, v2, v3, v4)] <= 0);
		model.addConstr(-zv[v3] + ze[e(v1, v3, v4)] + ze[e(v2, v3, v4)] - ze[e(v1, v2, v3, v4)] <= 0);
		model.addConstr(-zv[v4] + ze[e(v1, v3, v4)] + ze[e(v2, v3, v4)] - ze[e(v1, v2, v3, v4)] <= 0);

		model.addConstr(-zv[v1] - zv[v2] - zv[v3] - zv[v4] + 2 * ze[e(v1, v2)] + 2 * ze[e(v1, v3)] + 2 * ze[e(v1, v4)] + 2 * ze[e(v2, v3)] + 2 * ze[e(v2, v4)] + 2 * ze[e(v3, v4)] - 4 * ze[e(v1, v2, v3)] - 4 * ze[e(v1, v2, v4)] - 4 * ze[e(v1, v3, v4)] - 4 * ze[e(v2, v3, v4)] + 8 * ze[e(v1, v2, v3, v4)] + 1 == 1);
	}

	std::ofstream fin(output_file);

	std::cout << "runLP" << std::endl;

	for (int co = 0; co < 31; co++) {
		double corrupt = 0.00 + co * 0.01;
		std::cout << corrupt << std::endl;
		std::string temp = "corruption level: " + std::to_string(corrupt) + "\n";
		fin << temp;

		/* Record the tightness rate and recovery rate*/
		double tightness = 0;
		double partial_recovery = 0;
		double exact_recovery = 0;

		for (int num_sol = 0; num_sol < num_trials; num_sol++) {
			std::string path = filename;
			path += "code_" + std::to_string(num_bit) + "_" + std::to_string(corrupt) + "_" + std::to_string(num_sol) + ".txt";
			input_code->read_code(path, num_bit);
			GRBLinExpr obj = 1;

			for (int i = 0; i < num_bit; i++) {
				obj += (-2 * input_code->get_noisy_code()[i] + 1) * zv[i];
			}
			model.setObjective(obj, GRB_MINIMIZE);
			model.reset();
			auto begin = std::clock();
			model.optimize();
			auto end = std::clock();
			auto duration = double(end - begin) / CLOCKS_PER_SEC;
			// store results separated by space
			fin << " instance " << num_sol << ":";
			fin << duration << " ";
			double opt = obj.getValue();
			fin << opt << " ";

			double recovery_rate = 0;
			bool is_recovered = true;
			bool is_binary = true;
			for (int i = 0; i < num_bit; i++) {
				if (zv[i].get(GRB_DoubleAttr_X) > 0.0001 && zv[i].get(GRB_DoubleAttr_X) < 0.9999) {
					is_binary = false;
					is_recovered = false;
					break;
				}
			}
			if (is_binary) {
				for (int i = 0; i < num_bit; i++) {
					if (zv[i].get(GRB_DoubleAttr_X) < 0.0001) {
						recovery_rate += 1;
					}
					else {
						is_recovered = false;
					}
				}
			}
			else {
				for (int i = 0; i < num_bit; i++) {
					if (std::round(zv[i].get(GRB_DoubleAttr_X)) < 0.0001) {
						recovery_rate += 1;
					}
				}
			}

			fin << is_binary << " ";
			fin << recovery_rate / num_bit << " ";
			fin << is_recovered << " ";
			fin << "\n";

			tightness += (int)is_binary;
			exact_recovery += is_recovered;
			partial_recovery += recovery_rate;

			if (print_level >= 1) {
				std::cout << " clique: ";
				std::cout << " time: " << duration << " obj value: " << obj.getValue() << " gap: " << std::abs(obj.getValue() - opt) / std::abs(opt) << std::endl;
				if (print_level >= 2) {
					std::cout << " Coding Result: ";
					for (int i = 0; i < num_bit; i++) {
						std::cout << zv[i].get(GRB_DoubleAttr_X) << " ";
					}
					std::cout << std::endl;
				}
				std::cout << "   Is solution recovered: " << is_recovered << std::endl;

				if (print_level >= 3) {

					std::cout << " Received Code: ";
					for (int i = 0; i < num_bit; i++) {
						std::cout << input_code->get_noisy_code()[i] << " ";
					}
					std::cout << std::endl;

					std::cout << " Original Code: ";
					for (int i = 0; i < num_bit; i++) {
						std::cout << input_code->get_code()[i] << " ";
					}
					std::cout << std::endl;

					bool parity_check = true;
					for (const auto& row : input_code->get_check_matrix()) {
						int temp = 0;
						for (int i = 0; i < num_bit; i++) {
							temp += input_code->get_noisy_code()[i] * row[i];
						}
						if (temp % 2 != 0) {
							parity_check = false;
						}
					}

					std::cout << " Is noisy code satisfies the parity check: " << parity_check << std::endl;

					if (print_level >= 3) {
						std::cout << " Coding Result: ";
						for (int i = 0; i < num_bit; i++) {
							std::cout << zv[i].get(GRB_DoubleAttr_X) << " ";
						}
						std::cout << std::endl;
					}
				}
			}
		}

		std::cout << "tightness rate " << tightness / num_trials << " ";
		std::cout << "partial recovery rate " << partial_recovery / (num_trials * num_bit) << " ";
		std::cout << "exact recovery rate " << exact_recovery / num_trials << std::endl;
	}

	// clean up
	delete[] zv;
	delete[] ze;
}

void decoder::parLP(generator* input_code, int print_level, std::string output_file, std::string version, int num_trials) {

	int num_graph = input_code->get_check_matrix().size();
	int num_bit = input_code->get_check_matrix()[0].size();
	GRBEnv env = GRBEnv(true);
	GRBVar* zv = 0;
	env.set(GRB_IntParam_OutputFlag, 0);
	env.start();
	GRBModel model = GRBModel(env);
	model.set(GRB_IntParam_OutputFlag, 0);
	model.set(GRB_IntParam_LogToConsole, 0);

	zv = new GRBVar[num_bit];

	if (version == "LP") {
		for (int i = 0; i < num_bit; i++) {
			zv[i] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
		}
	}
	else {
		for (int i = 0; i < num_bit; i++) {
			zv[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		}
	}

	// parity lp
	for (const auto& row : input_code->get_check_matrix()) {
		std::vector<int> arr;
		for (int i = 0; i < row.size(); i++) {
			if (row[i] == 1) {
				arr.push_back(i);
			}
		}
		size_t n = arr.size();
		// Generate all subsets, particularly those with odd cardinality
		for (int subset = 1; subset < (1 << n); subset++) {
			int count = 0;
			std::vector<int> selectedIndices;

			for (int i = 0; i < n; i++) {
				if (subset & (1 << i)) {
					selectedIndices.push_back(i);
					count++;
				}
			}

			if (count % 2 != 0) { // Process only odd-sized subsets
				GRBLinExpr constr = 0;
				for (int i = 0; i < n; i++) {
					if (subset & (1 << i)) {
						constr += zv[arr[i]];
					}
					else {
						constr += -zv[arr[i]];
					}
				}
				int rhs = count - 1;
				model.addConstr(constr<=rhs);
			}
		}
	}

	std::ofstream fin(output_file);
	std::cout<<"Parity"<<std::endl;

	for (int co = 0; co < 21; co++) {
		double corrupt = 0.00 + co * 0.01;
		std::cout<<corrupt<<std::endl;
		std::string temp = "corruption level: " + std::to_string(corrupt) + "\n";
		fin << temp;

		/* Record the tightness rate and recovery rate*/
		double tightness = 0;
		double partial_recovery = 0;
		double exact_recovery = 0;

		for (int num_sol = 0; num_sol < num_trials; num_sol++) {
			std::string path = filename;
			path += "code_" + std::to_string(num_bit) + "_" + std::to_string(corrupt) + "_" + std::to_string(num_sol) + ".txt";
			input_code->read_code(path, num_bit);
			GRBLinExpr obj = 1;

			for (int i = 0; i < num_bit; i++) {
				obj += (-2 * input_code->get_noisy_code()[i] + 1) * zv[i];
			}
			model.setObjective(obj, GRB_MINIMIZE);

			auto begin = std::clock();
			model.optimize();
			auto end = std::clock();
			auto duration = double(end - begin) / CLOCKS_PER_SEC;

			fin << " instance " << num_sol << ":";
			fin << duration << " ";
			double opt = obj.getValue();
			fin << opt << " ";

			double recovery_rate = 0;
			bool is_recovered = true;
			bool is_binary = true;
			for (int i = 0; i < num_bit; i++) {
				if (zv[i].get(GRB_DoubleAttr_X) > 0.0001 && zv[i].get(GRB_DoubleAttr_X) < 0.9999) {
					is_binary = false;
					is_recovered = false;
					break;
				}
			}
			if (is_binary) {
				for (int i = 0; i < num_bit; i++) {
					if (zv[i].get(GRB_DoubleAttr_X) < 0.0001) {
						recovery_rate += 1;
					}
					else {
						is_recovered = false;
					}
				}
			}
			else {
				for (int i = 0; i < num_bit; i++) {
					if (std::round(zv[i].get(GRB_DoubleAttr_X)) < 0.0001) {
						recovery_rate += 1;
					}
				}
			}

			fin << is_binary << " ";
			fin << recovery_rate / num_bit << " ";
			fin << is_recovered << " ";
			fin << "\n";

			tightness += (int)is_binary;
			exact_recovery += is_recovered;
			partial_recovery += recovery_rate;

			if (print_level >= 1) {
				std::cout << " parity IP: ";
				std::cout << " time: " << duration << " obj value: " << obj.getValue() << " gap: " << std::abs(obj.getValue() - opt) / std::abs(opt) << std::endl;
				if (print_level >= 2) {
					std::cout << " Coding Result: ";
					for (int i = 0; i < num_bit; i++) {
						std::cout << zv[i].get(GRB_DoubleAttr_X) << " ";
					}
					std::cout << std::endl;
				}
				std::cout << "   Is solution recovered: " << is_recovered << std::endl;

				if (print_level >= 2) {

					std::cout << " Received Code: ";
					for (int i = 0; i < num_bit; i++) {
						std::cout << input_code->get_noisy_code()[i] << " ";
					}
					std::cout << std::endl;

					std::cout << " Original Code: ";
					for (int i = 0; i < num_bit; i++) {
						std::cout << input_code->get_code()[i] << " ";
					}
					std::cout << std::endl;

					bool parity_check = true;
					for (const auto& row : input_code->get_check_matrix()) {
						int temp = 0;
						for (int i = 0; i < num_bit; i++) {
							temp += input_code->get_noisy_code()[i] * row[i];
						}
						if (temp % 2 != 0) {
							parity_check = false;
						}
					}

					std::cout << " Is noisy code satisfies the parity check: " << parity_check << std::endl;

					if (print_level >= 3) {
						std::cout << " Coding Result: ";
						for (int i = 0; i < num_bit; i++) {
							std::cout << zv[i].get(GRB_DoubleAttr_X) << " ";
						}
						std::cout << std::endl;
					}
				}
			}
		}

		std::cout << "tightness rate " << tightness / num_trials << " ";
		std::cout << "partial recovery rate " << partial_recovery / (num_trials * num_bit) << " ";
		std::cout << "exact recovery rate " << exact_recovery / num_trials << std::endl;
	}

	// clean up
	delete[] zv;
}

void decoder::cliqueLP(generator* input_code, int print_level, std::string output_file, std::string version, int num_trials) {

	int num_graph = input_code->get_check_matrix().size();
	int num_bit = input_code->get_check_matrix()[0].size();

	GRBEnv env = GRBEnv(true);
	GRBVar* zv = 0;
	GRBVar* ze = 0;
	env.set(GRB_IntParam_OutputFlag, 0);
	env.start();
	GRBModel model = GRBModel(env);
	model.set(GRB_IntParam_OutputFlag, 0);
	model.set(GRB_IntParam_LogToConsole, 0);

	zv = new GRBVar[num_bit];
	ze = new GRBVar[mapping.size()];


	if (version == "LP") {
		for (int i = 0; i < num_bit; i++) {
			zv[i] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
		}
		for (int i = 0; i < mapping.size(); i++) {
			ze[i] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
		}
	}
	else {
		for (int i = 0; i < num_bit; i++) {
			zv[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		}
		for (int i = 0; i < mapping.size(); i++) {
			ze[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		}
	}



	// clique lp
	for (const auto& row : input_code->get_check_matrix()) {
		std::vector<int> arr;
		for (int i = 0; i < row.size(); i++) {
			if (row[i] == 1) {
				arr.push_back(i);
			}
		}
		GRBLinExpr parity_hyperplane;

		size_t n = arr.size();
		// Generate all subsets, particularly those with odd cardinality
		for (int subset = 1; subset < (1 << n); subset++) {
			std::vector<int> selectedIndices;
			std::vector<int> remainingIndices;
			std::vector<int> parity_hyperplane_term;
			// Check each bit position
			for (int i = 0; i < n; i++) {
				if (subset & (1 << i)) {
					selectedIndices.push_back(i);
					parity_hyperplane_term.push_back(arr[i]);
				}
				else {
					remainingIndices.push_back(i);
				}
			}

			/*Add parity constraints*/
			if (parity_hyperplane_term.size() == 1) {
				parity_hyperplane += std::pow(2, (parity_hyperplane_term.size() - 1)) * (parity_hyperplane_term.size() % 2 == 1 ? -1 : 1) * zv[parity_hyperplane_term[0]];
			}
			else {
				parity_hyperplane += std::pow(2, (parity_hyperplane_term.size() - 1)) * (parity_hyperplane_term.size() % 2 == 1 ? -1 : 1) * ze[edge(parity_hyperplane_term)];
			}



			std::vector<std::pair<std::vector<int>, int>> terms{ {{}, 1} }; // Start with an empty product and a coefficient of +1
			for (int idx : selectedIndices) {
				std::vector<std::pair<std::vector<int>, int>> newTerms;
				for (auto& term : terms) {
					// Add the current term (+1)
					newTerms.push_back({ term.first, term.second });
					// Create a new term with the current value (-zi), flipping the coefficient
					std::vector<int> newValues = term.first;
					newValues.push_back(arr[idx]);
					newTerms.push_back({ newValues, -term.second });
				}
				terms = newTerms;
			}

			// Multiply by remaining values
			if (!remainingIndices.empty()) {
				for (auto& term : terms) {
					for (int idx : remainingIndices) {
						term.first.push_back(arr[idx]);
					}
				}
			}

			GRBLinExpr constr;

			for (auto& term : terms) {

				if (term.first.size() == 1) {
					constr += term.second * zv[term.first[0]];
				}
				else if (term.first.empty()) {
					constr += 1;
				}
				else {
					constr += term.second * ze[edge(term.first)];
				}
			}

			model.addConstr(constr >= 0);
		}

		model.addConstr(parity_hyperplane == 0);
	}

	std::ofstream fin(output_file);

	std::cout<< "Clique" << std::endl;

	for (int co = 0; co < 21; co++) {
		double corrupt = 0.00 + co * 0.01;
		std::cout << corrupt << std::endl;
		std::string temp = "corruption level: " + std::to_string(corrupt) + "\n";
		fin << temp;

		/* Record the tightness rate and recovery rate*/
		double tightness = 0;
		double partial_recovery = 0;
		double exact_recovery = 0;

		for (int num_sol = 0; num_sol < num_trials; num_sol++) {
			std::string path = filename;
			path += "code_" + std::to_string(num_bit) + "_" + std::to_string(corrupt) + "_" + std::to_string(num_sol) + ".txt";
			input_code->read_code(path, num_bit);
			GRBLinExpr obj = 1;

			for (int i = 0; i < num_bit; i++) {
				obj += (-2 * input_code->get_noisy_code()[i] + 1) * zv[i];
			}
			model.setObjective(obj, GRB_MINIMIZE);
			model.reset();
			auto begin = std::clock();
			model.optimize();
			auto end = std::clock();
			auto duration = double(end - begin) / CLOCKS_PER_SEC;
			// store results separated by space
			fin << " instance " << num_sol << ":";
			fin << duration << " ";
			double opt = obj.getValue();
			fin << opt << " ";

			double recovery_rate = 0;
			bool is_recovered = true;
			bool is_binary = true;
			for (int i = 0; i < num_bit; i++) {
				if (zv[i].get(GRB_DoubleAttr_X) > 0.0001 && zv[i].get(GRB_DoubleAttr_X) < 0.9999) {
					is_binary = false;
					is_recovered = false;
					break;
				}
			}
			if (is_binary) {
				for (int i = 0; i < num_bit; i++) {
					if (zv[i].get(GRB_DoubleAttr_X) < 0.0001) {
						recovery_rate += 1;
					}
					else {
						is_recovered = false;
					}
				}
			}
			else {
				for (int i = 0; i < num_bit; i++) {
					if (std::round(zv[i].get(GRB_DoubleAttr_X)) < 0.0001) {
						recovery_rate += 1;
					}
				}
			}

			fin << is_binary << " ";
			fin << recovery_rate / num_bit << " ";
			fin << is_recovered << " ";
			fin << "\n";

			tightness += (int)is_binary;
			exact_recovery += is_recovered;
			partial_recovery += recovery_rate;

			if (print_level >= 1) {
				std::cout << " clique: ";
				std::cout << " time: " << duration << " obj value: " << obj.getValue() << " gap: " << std::abs(obj.getValue() - opt) / std::abs(opt) << std::endl;
				if (print_level >= 2) {
					std::cout << " Coding Result: ";
					for (int i = 0; i < num_bit; i++) {
						std::cout << zv[i].get(GRB_DoubleAttr_X) << " ";
					}
					std::cout << std::endl;
				}
				std::cout << "   Is solution recovered: " << is_recovered << std::endl;

				if (print_level >= 3) {

					std::cout << " Received Code: ";
					for (int i = 0; i < num_bit; i++) {
						std::cout << input_code->get_noisy_code()[i] << " ";
					}
					std::cout << std::endl;

					std::cout << " Original Code: ";
					for (int i = 0; i < num_bit; i++) {
						std::cout << input_code->get_code()[i] << " ";
					}
					std::cout << std::endl;

					bool parity_check = true;
					for (const auto& row : input_code->get_check_matrix()) {
						int temp = 0;
						for (int i = 0; i < num_bit; i++) {
							temp += input_code->get_noisy_code()[i] * row[i];
						}
						if (temp % 2 != 0) {
							parity_check = false;
						}
					}

					std::cout << " Is noisy code satisfies the parity check: " << parity_check << std::endl;

					if (print_level >= 3) {
						std::cout << " Coding Result: ";
						for (int i = 0; i < num_bit; i++) {
							std::cout << zv[i].get(GRB_DoubleAttr_X) << " ";
						}
						std::cout << std::endl;
					}
				}
			}
		}

		std::cout << "tightness rate " << tightness / num_trials << " ";
		std::cout << "partial recovery rate " << partial_recovery / (num_trials * num_bit) << " ";
		std::cout << "exact recovery rate " << exact_recovery / num_trials << std::endl;
	}

	// clean up
	delete[] zv;
	delete[] ze;
}

void decoder::McliqueLP(generator* input_code, int print_level, std::string output_file, std::string version, int cycle_size, int num_trials) {
	int num_graph = input_code->get_check_matrix().size();
	int num_bit = input_code->get_check_matrix()[0].size();

	GRBEnv env = GRBEnv(true);
	GRBVar* zv = 0;
	GRBVar* ze = 0;
	env.set(GRB_IntParam_OutputFlag, 0);
	env.start();
	GRBModel model = GRBModel(env);
	model.set(GRB_IntParam_OutputFlag, 0);
	model.set(GRB_IntParam_LogToConsole, 0);

	zv = new GRBVar[num_bit];
	ze = new GRBVar[mapping.size()];


	if (version == "LP") {
		for (int i = 0; i < num_bit; i++) {
			zv[i] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
		}
		for (int i = 0; i < mapping.size(); i++) {
			ze[i] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
		}
	}
	else {
		for (int i = 0; i < num_bit; i++) {
			zv[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		}
		for (int i = 0; i < mapping.size(); i++) {
			ze[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		}
	}



	// clique lp
	for (const auto& row : input_code->get_check_matrix()) {
		std::vector<int> arr;
		for (int i = 0; i < row.size(); i++) {
			if (row[i] == 1) {
				arr.push_back(i);
			}
		}
		GRBLinExpr parity_hyperplane;

		size_t n = arr.size();
		// Generate all subsets, particularly those with odd cardinality
		for (int subset = 1; subset < (1 << n); subset++) {
			std::vector<int> selectedIndices;
			std::vector<int> remainingIndices;
			std::vector<int> parity_hyperplane_term;
			// Check each bit position
			for (int i = 0; i < n; i++) {
				if (subset & (1 << i)) {
					selectedIndices.push_back(i);
					parity_hyperplane_term.push_back(arr[i]);
				}
				else {
					remainingIndices.push_back(i);
				}
			}

			/*Add parity constraints*/
			if (parity_hyperplane_term.size() == 1) {
				parity_hyperplane += std::pow(2, (parity_hyperplane_term.size() - 1)) * (parity_hyperplane_term.size() % 2 == 1 ? -1 : 1) * zv[parity_hyperplane_term[0]];
			}
			else {
				parity_hyperplane += std::pow(2, (parity_hyperplane_term.size() - 1)) * (parity_hyperplane_term.size() % 2 == 1 ? -1 : 1) * ze[edge(parity_hyperplane_term)];
			}



			std::vector<std::pair<std::vector<int>, int>> terms{ {{}, 1} }; // Start with an empty product and a coefficient of +1
			for (int idx : selectedIndices) {
				std::vector<std::pair<std::vector<int>, int>> newTerms;
				for (auto& term : terms) {
					// Add the current term (+1)
					newTerms.push_back({ term.first, term.second });
					// Create a new term with the current value (-zi), flipping the coefficient
					std::vector<int> newValues = term.first;
					newValues.push_back(arr[idx]);
					newTerms.push_back({ newValues, -term.second });
				}
				terms = newTerms;
			}

			// Multiply by remaining values
			if (!remainingIndices.empty()) {
				for (auto& term : terms) {
					for (int idx : remainingIndices) {
						term.first.push_back(arr[idx]);
					}
				}
			}

			GRBLinExpr constr;

			for (auto& term : terms) {

				if (term.first.size() == 1) {
					constr += term.second * zv[term.first[0]];
				}
				else if (term.first.empty()) {
					constr += 1;
				}
				else {
					constr += term.second * ze[edge(term.first)];
				}
			}

			model.addConstr(constr >= 0);
		}

		model.addConstr(parity_hyperplane == 0);
	}

	// multi clique inequality

	std::vector<GRBLinExpr> cycle_ineq1;
	std::vector<GRBLinExpr> lifted_cycle_ineq1;
	std::vector<GRBLinExpr> lifted_cycle_ineq2;

	int structure_count = 0;

	for (int size = 3; size <= cycle_size; size++) {
		// size 3 cycles
		std::vector<std::vector<int>> combinations;
		std::vector<int> currentCombination;

		// Generate subset of cliques with size 'size'
		generateCombinations(currentCombination, combinations, num_graph, size, 0);


		std::vector<std::vector<int>> nonzerocol;
		for (const auto& row : input_code->get_check_matrix()) {
			std::vector<int> subgraph1;
			for (int col = 0; col < row.size(); col++) {
				if (row[col] != 0) {
					subgraph1.push_back(col);
				}
			}
			nonzerocol.push_back(subgraph1);
		}

		for (const auto& row : combinations) {
			bool is_structure = false;
			std::vector<std::vector<int>> sub_graph;
			for (int i = 0; i < size; i++) {
				sub_graph.push_back(nonzerocol[row[i]]);
			}

			std::vector<int> cycle;

			if (std::find_first_of(sub_graph[0].begin(), sub_graph[0].end(),
				sub_graph[1].begin(), sub_graph[1].end()) != sub_graph[0].end()) {

				int commonElement = -10;
				std::vector<std::vector<int>> modifiedVectors;

				// lifted cycle ineq
				if (findCommonElementAndExclude(sub_graph, commonElement, modifiedVectors)) {
					is_structure = Check_Structure(modifiedVectors, cycle);
					if (is_structure) {
						structure_count++;

						std::vector<std::vector<int>> edges;

						// Create the full cycle of edges
						for (size_t i = 0; i < cycle.size(); i++) {
							edges.push_back({ cycle[i], cycle[(i + 1) % cycle.size()] });
						}

						int n = edges.size();
						// Process each subset
						for (int subset = 1; subset < (1 << n); subset++) {

							int count = 0;
							for (int i = 0; i < n; i++) {
								if (subset & (1 << i)) {
									count++;
								}
							}

							GRBLinExpr cycle_ineq1 = 0;
							GRBLinExpr cycle_ineq2 = 0;

							if (count % 2 == 1) { // Check if the subset has an odd number of edges
								std::vector<std::vector<int>> S1, S2;
								std::unordered_set<int> nodesInS1, nodesInS2;

								for (int i = 0; i < n; ++i) {
									if (subset & (1 << i)) { // Check if the ith edge is in the subset
										S1.push_back(edges[i]);
										nodesInS1.insert(edges[i][0]);
										nodesInS1.insert(edges[i][1]);
									}
									else {
										S2.push_back(edges[i]);
										nodesInS2.insert(edges[i][0]);
										nodesInS2.insert(edges[i][1]);
									}
								}

								// add terms to the inequality

								for (auto& item : S1) {
									// add common element to the item
									std::vector<int> new_item = item;
									new_item.push_back(commonElement);
									cycle_ineq1 += -ze[edge(new_item)];
									cycle_ineq2 += -ze[edge(item)];
									cycle_ineq2 += ze[edge(new_item)];
								}

								for (auto& item : S2) {
									std::vector<int> new_item = item;
									new_item.push_back(commonElement);
									cycle_ineq1 += ze[edge(new_item)];
									cycle_ineq2 += ze[edge(item)];
									cycle_ineq2 += -ze[edge(new_item)];
								}

								// Print nodes based on their presence in S1 and S2
								for (int node : cycle) {
									if (nodesInS2.find(node) == nodesInS2.end()) {
										cycle_ineq1 += ze[edge(std::vector<int>{node, commonElement})];
										cycle_ineq2 += zv[node];
										cycle_ineq2 += -ze[edge(std::vector<int>{node, commonElement})];
									}
									else if (nodesInS1.find(node) == nodesInS1.end()) {
										cycle_ineq1 += -ze[edge(std::vector<int>{node, commonElement})];
										cycle_ineq2 += -zv[node];
										cycle_ineq2 += ze[edge(std::vector<int>{node, commonElement})];
									}
								}
								cycle_ineq1 += -std::floor(S1.size() / 2)*zv[commonElement];
								cycle_ineq2 += -std::floor(S1.size() / 2)*(1- zv[commonElement]);
							}

							lifted_cycle_ineq1.push_back(cycle_ineq1);
							lifted_cycle_ineq2.push_back(cycle_ineq2);
							//model.addConstr(cycle_ineq1 <= 0);
							//model.addConstr(cycle_ineq2 <= 0);
						}

					}
				}
				// odd cycle
				else {
					is_structure = Check_Structure(sub_graph, cycle);

					if (is_structure) {
						structure_count++;

						std::vector<std::vector<int>> edges;

						// Create the full cycle of edges
						for (size_t i = 0; i < cycle.size(); i++) {
							edges.push_back({ cycle[i], cycle[(i + 1) % cycle.size()] });
						}

						int n = edges.size();
						// Process each subset
						for (int subset = 1; subset < (1 << n); subset++) {

							int count = 0;
							for (int i = 0; i < n; i++) {
								if (subset & (1 << i)) {
									count++;
								}
							}

							GRBLinExpr cycle_ineq = 0;

							if (count % 2 == 1) { // Check if the subset has an odd number of edges
								std::vector<std::vector<int>> S1, S2;
								std::unordered_set<int> nodesInS1, nodesInS2;

								for (int i = 0; i < n; ++i) {
									if (subset & (1 << i)) { // Check if the ith edge is in the subset
										S1.push_back(edges[i]);
										nodesInS1.insert(edges[i][0]);
										nodesInS1.insert(edges[i][1]);
									}
									else {
										S2.push_back(edges[i]);
										nodesInS2.insert(edges[i][0]);
										nodesInS2.insert(edges[i][1]);
									}
								}

								for (auto& item : S1) {
									cycle_ineq += -ze[edge(item)];
								}

								for (auto& item : S2) {
									cycle_ineq += ze[edge(item)];
								}

								// Print nodes based on their presence in S1 and S2
								for (int node : cycle) {
									if (nodesInS2.find(node) == nodesInS2.end()) {
										cycle_ineq += zv[node];
									}
									else if (nodesInS1.find(node) == nodesInS1.end()) {
										cycle_ineq += -zv[node];
									}
								}
								cycle_ineq += -std::floor(S1.size() / 2);
							}

							cycle_ineq1.push_back(cycle_ineq);
							//model.addConstr(cycle_ineq <= 0);
						}

					}
				}
			}
			else {
				continue;
			}
		}
	}

	std::cout<<"structure count "<<structure_count<<std::endl;

	std::ofstream fin(output_file);

	std::cout << "Multi-Clique" << std::endl;

	for (int co = 0; co < 21; co++) {
		double corrupt = 0.00 + co * 0.01;
		std::cout << corrupt << std::endl;
		std::string temp = "corruption level: " + std::to_string(corrupt) + "\n";
		fin << temp;

		/* Record the tightness rate and recovery rate*/
		double tightness = 0;
		double partial_recovery = 0;
		double exact_recovery = 0;

		for (int num_sol = 0; num_sol < num_trials; num_sol++) {
			std::string path = filename;
			path += "code_" + std::to_string(num_bit) + "_" + std::to_string(corrupt) + "_" + std::to_string(num_sol) + ".txt";
			input_code->read_code(path, num_bit);
			GRBLinExpr obj = 1;

			for (int i = 0; i < num_bit; i++) {
				obj += (-2 * input_code->get_noisy_code()[i] + 1) * zv[i];
			}
			model.setObjective(obj, GRB_MINIMIZE);

			auto begin = std::clock();
			model.optimize();
			auto end = std::clock();
			auto duration = double(end - begin) / CLOCKS_PER_SEC;


			// store results separated by space
			double opt = 0.0;
			double recovery_rate = 0;
			bool is_recovered = true;
			bool is_binary = true;

			for (int i = 0; i < num_bit; i++) {
				if (zv[i].get(GRB_DoubleAttr_X) > 0.0001 && zv[i].get(GRB_DoubleAttr_X) < 0.9999) {
					is_binary = false;
					break;
				}
			}
			if (is_binary) {
				for (int i = 0; i < num_bit; i++) {
					if (zv[i].get(GRB_DoubleAttr_X) >= 0.5) {
						is_recovered = false;
					}
					if (zv[i].get(GRB_DoubleAttr_X) < 0.0001) {
						recovery_rate += 1;
					}
				}
			}

			// if non-binary solution, check the lifted cycle inequalities

			std::vector<GRBConstr> lifted_cycle_constr;
			// resize
			lifted_cycle_constr.resize(lifted_cycle_ineq1.size()+ lifted_cycle_ineq2.size()+ cycle_ineq1.size());
			int num_new_cons = 0;


			if (!is_binary) {
				// add lifted cycle inequalities and resolve
				for (int i = 0; i < lifted_cycle_ineq1.size(); i++) {
					if (lifted_cycle_ineq1[i].getValue() > 1e-3) {
						lifted_cycle_constr[num_new_cons] = model.addConstr(lifted_cycle_ineq1[i] <= 0);
						num_new_cons++;
					}
					if (lifted_cycle_ineq2[i].getValue() > 1e-3) {
						lifted_cycle_constr[num_new_cons] = model.addConstr(lifted_cycle_ineq2[i] <= 0);
						num_new_cons++;
					}
					//lifted_cycle_constr[i] = model.addConstr(lifted_cycle_ineq1[i] <= 0);
					//lifted_cycle_constr[i+ lifted_cycle_ineq1.size()] = model.addConstr(lifted_cycle_ineq2[i] <= 0);
				}
				for (int i = 0; i < cycle_ineq1.size(); i++) {
					if (cycle_ineq1[i].getValue() > 1e-3) {
						lifted_cycle_constr[num_new_cons] = model.addConstr(cycle_ineq1[i] <= 0);
						num_new_cons++;
					}
					//lifted_cycle_constr[i+2*lifted_cycle_ineq1.size()] = model.addConstr(cycle_ineq1[i] <= 0);
				}
				
				auto begin = std::clock();
				model.optimize();
				auto end = std::clock();
				auto duration = double(end - begin) / CLOCKS_PER_SEC;

				is_recovered = true;
				is_binary = true;
				recovery_rate = 0;
				for (int i = 0; i < num_bit; i++) {
					if (zv[i].get(GRB_DoubleAttr_X) > 0.0001 && zv[i].get(GRB_DoubleAttr_X) < 0.9999) {
						is_binary = false;
						is_recovered = false;
						break;
					}
				}
				if (is_binary) {
					for (int i = 0; i < num_bit; i++) {
						if (zv[i].get(GRB_DoubleAttr_X) >= 0.5) {
							is_recovered = false;
						}
						if (zv[i].get(GRB_DoubleAttr_X) < 0.0001) {
							recovery_rate += 1;
						}
					}
				}
				else {
					for (int i = 0; i < num_bit; i++) {
						if (std::round(zv[i].get(GRB_DoubleAttr_X)) < 0.0001) {
							recovery_rate += 1;
						}
					}
				}

				// store results separated by space
				fin << " instance " << num_sol << ":";
				fin << duration << " ";
				opt = obj.getValue();
				fin << opt << " ";

				fin << is_binary << " ";
				fin << recovery_rate / num_bit << " ";
				fin << is_recovered << " ";
				fin << "\n";

				tightness += (int)is_binary;
				exact_recovery += is_recovered;
				partial_recovery += recovery_rate;

				// remove constraints
				//for (int i = 0; i < lifted_cycle_constr.size(); i++) {
				//	model.remove(lifted_cycle_constr[i]);
				//}
				for (int i = 0; i < num_new_cons; i++) {
					model.remove(lifted_cycle_constr[i]);
				}

				model.update();
			}
			else {
				// store results separated by space
				fin << " instance " << num_sol << ":";
				fin << duration << " ";
				opt = obj.getValue();
				fin << opt << " ";

				fin << is_binary << " ";
				fin << recovery_rate / num_bit << " ";
				fin << is_recovered << " ";
				fin << "\n";

				tightness += (int)is_binary;
				exact_recovery += is_recovered;
				partial_recovery += recovery_rate;
			}

			lifted_cycle_constr.clear();

			if (print_level >= 1) {
				std::cout << " multi-cliques: ";
				std::cout << " time: " << duration << " obj value: " << obj.getValue() << " gap: " << std::abs(obj.getValue() - opt) / std::abs(opt) << std::endl;
				if (print_level >= 2) {
					std::cout << " Coding Result: ";
					for (int i = 0; i < num_bit; i++) {
						std::cout << zv[i].get(GRB_DoubleAttr_X) << " ";
					}
					std::cout << std::endl;
				}
				std::cout << "   Is solution recovered: " << is_recovered << std::endl;

				if (print_level >= 3) {

					std::cout << " Received Code: ";
					for (int i = 0; i < num_bit; i++) {
						std::cout << input_code->get_noisy_code()[i] << " ";
					}
					std::cout << std::endl;

					std::cout << " Original Code: ";
					for (int i = 0; i < num_bit; i++) {
						std::cout << input_code->get_code()[i] << " ";
					}
					std::cout << std::endl;

					bool parity_check = true;
					for (const auto& row : input_code->get_check_matrix()) {
						int temp = 0;
						for (int i = 0; i < num_bit; i++) {
							temp += input_code->get_noisy_code()[i] * row[i];
						}
						if (temp % 2 != 0) {
							parity_check = false;
						}
					}

					std::cout << " Is noisy code satisfies the parity check: " << parity_check << std::endl;

					if (print_level >= 3) {
						std::cout << " Coding Result: ";
						for (int i = 0; i < num_bit; i++) {
							std::cout << zv[i].get(GRB_DoubleAttr_X) << " ";
						}
						std::cout << std::endl;
					}
				}
			}
		}

		std::cout << "tightness rate " << tightness / num_trials << " ";
		std::cout << "partial recovery rate " << partial_recovery / (num_trials * num_bit) << " ";
		std::cout << "exact recovery rate " << exact_recovery / num_trials << std::endl;
	}

	// clean up
	delete[] zv;
	delete[] ze;

}