#pragma once
#include "Generator.h"


class decoder {
public:
	decoder(generator* input_code, std::string filename);
	void parLP(generator* input_code, int print_level, std::string output_file, std::string version, int num_trials);
	void cliqueLP(generator* input_code, int print_level, std::string output_file, std::string version, int num_trials);
	void stdLP(generator* input_code, int print_level, std::string output_file, std::string version, int num_trials);
	void flLP(generator* input_code, int print_level, std::string output_file, std::string version, int num_trials);
	void runLP(generator* input_code, int print_level, std::string output_file, std::string version, int num_trials);
	void McliqueLP(generator* input_code, int print_level, std::string output_file, std::string version, int cycle_size, int num_trials);
private:
	std::string filename;
	std::unordered_map<std::string, int> mapping;
	int e(int i = -1, int j = -1, int k = -1, int l = -1);
	int edge(std::vector<int> arr);
	bool Check_Structure(std::vector<std::vector<int>> Three_HyperGraphs, std::vector<int>& cycle);
	bool findCycle(std::vector<std::vector<int>>& sets, std::vector<int>& visited_sets, std::vector<int>& cycle, int startSet, int currentSet);
	void generateCombinations(std::vector<int>& currentCombination, std::vector<std::vector<int>>& combinations, int n, int k, int start);
	bool findCommonElementAndExclude(const std::vector<std::vector<int>>& vectors, int& commonElement, std::vector<std::vector<int>>& modifiedVectors);
};