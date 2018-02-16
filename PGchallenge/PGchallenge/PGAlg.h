#pragma once
#include <windows.h>
#include <stdlib.h>
#include <atlstr.h>
#include "CTree.h"
#include <cmath>
#include <ctime>
#include <vector>
#include <fstream>
#include <iomanip>
#include <functional>  // sorting
#include <iostream>  // debug

static const int i_NUMBER_OF_VARS = 2;
static const int i_SAMPLE_DATA_SIZE = 1000;
static const double d_VARAIBLE_MAX_VAL = 100;
static const double d_VARIABLE_MIN_VAL = -100;
static const std::string s_FUNCTION_FILE_NAME = "model_function.txt";
static const std::string s_DATA_FILE_NAME = "model_data.txt";
static const char c_DELIMITER = ';';
static const double d_EPSILON_ACCURACY = 1e-10;
static const int i_WRITE_TO_FILE_PRECISION = 10;
static const int i_TREE_MAX_HEIGHT = 3;
static const int i_INITIAL_POPULATION_SIZE = 100;
static const int i_NUM_OF_ITERATIONS = 100;
static const int i_MIX_IN_BESTS_INTERVAL = i_NUM_OF_ITERATIONS / 4;
static const int i_ACCEPTABLE_ITERATIONS_WITHOUT_IMPROVEMENT = i_NUM_OF_ITERATIONS / 20;
static const int i_LINEAR_REGRESSION_SAMPLE = i_NUM_OF_ITERATIONS / 10;
static const double d_CROSSOVER_PROBABILITY = 0.6;
static const double d_MUTATION_PROBABILITY = 0.4;
static const double d_INCREASED_CROSSOVER_PROB = 0.7;
static const double d_INCREASED_MUTATION_PROB = 0.9;
static const int i_HOW_MANY_TO_MIXIN = i_NUM_OF_ITERATIONS / 4;
static const int i_HOW_MANY_TO_REPLACE_WTH_NEW = i_INITIAL_POPULATION_SIZE / 2;
static const int i_ACCEPTABLE_NUM_OF_NODES = std::pow(2, i_TREE_MAX_HEIGHT + 1) * 2;
static const int i_POPULATION_MOD_INTERVAL = i_LINEAR_REGRESSION_SAMPLE / 2;
static const int i_ITERATIONS_WITH_INCREASED_POPULATION = i_NUM_OF_ITERATIONS / 20;
static const double d_INCREASE_POPULATION_FACTOR = 1.25;
static const int i_BEST_TREE_FAMILY_SIZE = i_NUM_OF_ITERATIONS / 4;

class CPGAlg
{
public:

	CPGAlg();
	~CPGAlg();

	bool bGenerateRndTreeFile(int iTreeMaxHeight);
	bool bGenereateFromStringTreeFile(std::string sFormula);
	bool bInitialize(CString sTest);
	void vRunIteration();
	CString sGetCurrentBestTree();
	double dGetBestDistance();

	void vPrintDataArray();  // debug, testing

private:

	bool b_initialized;
	std::vector<std::string> vs_split_string(std::string sToSplit, char cDelimiter);
	void v_initialize_rnd_seed();
	void v_delete_population();
	void v_delete_parents();
	void v_delete_best_tree();
	void v_evaluate_single_tree(CTree * c_tree);
	void v_evaluate_trees();
	void v_selection();
	void v_parents_crossover();
	void v_mutations();
	void v_mixin_bests();  // replaces N worst trees in population with bests
	void v_check_modifications();
	void v_check_normalize_params();
	void v_increase_probabilities();
	void v_normalize_probabilities();
	bool b_nodes_num_optimization_needed();
	void v_replace_worsts_with_new();
	void v_increase_population_size();
	void v_normalize_population_size();
	void v_mixin_mutated_best_tree_family();

	// member fields
	double d_mut_prob;
	double d_cross_prob;
	std::vector<std::vector<std::double_t>> vd_data;
	std::string s_current_best_tree;
	std::vector<CTree *> vc_population;
	std::vector<CTree *> vc_parents;
	CTree * c_current_best;
	double d_current_min_distance;
	std::vector<CTree*> vc_best_trees;
	int i_iteration_num;
	int i_iterations_without_improve;
	int i_modifications_iteration;
	int i_average_num_of_nodes;
	double d_average_distance;
	int i_population_mod_iteration;
	int i_iteration_pop_incr;
	bool b_population_increased;

	// linear regression (number of nodes as x and average distance as y), decreasing number of nodes decreases evaluaton time
	std::vector<int> vi_averages_nodes_num;
	std::vector<double> vd_averages_distances;

	static bool b_srand_was_called;
};// class CPGAlg


