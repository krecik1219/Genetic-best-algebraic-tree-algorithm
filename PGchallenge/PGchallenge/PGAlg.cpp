#include "PGAlg.h"

CPGAlg::CPGAlg()
{
	v_initialize_rnd_seed();
	c_current_best = NULL;
	d_current_min_distance = INFINITY;
	b_initialized = false;
	i_iteration_num = 0;
	d_mut_prob = d_MUTATION_PROBABILITY;
	d_cross_prob = d_CROSSOVER_PROBABILITY;
	i_iterations_without_improve = 0;
	i_modifications_iteration = 0;
	i_population_mod_iteration = 0;
	i_iteration_pop_incr = 0;
}// CPGAlg::CPGAlg()

CPGAlg::~CPGAlg()
{
	v_delete_best_tree();
	v_delete_population();
	v_delete_parents();
	vc_parents.clear();
	vc_population.clear();
	int i_best_trees_num = vc_best_trees.size();
	for (int ii = 0; ii < i_best_trees_num; ii++)
	{
		delete vc_best_trees[ii];
	}// for (int ii = 0; ii < i_best_trees_num; ii++)
	vc_best_trees.clear();
}// CPGAlg::~CPGAlg()

bool CPGAlg::bGenerateRndTreeFile(int iTreeMaxHeight)
{
	//srand(time(NULL));
	CTree * c_tree = new CTree(iTreeMaxHeight);
	double ** pd_model_data = new double *[i_NUMBER_OF_VARS + 1];
	for (int ii = 0; ii < i_NUMBER_OF_VARS + 1; ii++)
	{
		pd_model_data[ii] = new double[i_SAMPLE_DATA_SIZE];
	}// for (int ii = 0; ii < i_NUMBER_OF_VARS; ii++)

	double d_rnd_value;
	std::vector<double> d_eval_data(i_NUMBER_OF_VARS);
	for (int ii = 0; ii < i_SAMPLE_DATA_SIZE; ii++)
	{
		for (int ij = 0; ij < i_NUMBER_OF_VARS; ij++)
		{
			d_rnd_value = (double)rand() / RAND_MAX;
			d_rnd_value = d_VARIABLE_MIN_VAL + d_rnd_value * (d_VARAIBLE_MAX_VAL - d_VARIABLE_MIN_VAL);
			//d_rnd_value = round(d_rnd_value * 100000) / 1000.0;  // to avoid some problembs but it doesnt really work
			pd_model_data[ij][ii] = d_rnd_value; // variable
			d_eval_data[ij] = d_rnd_value;
		}// for (int ij = 0; ij < i_NUMBER_OF_VARS; ij++)
		pd_model_data[i_NUMBER_OF_VARS][ii] = c_tree->dEvaluateFormula(d_eval_data);
	}// for (int ii = 0; ii < i_SAMPLE_DATA_SIZE; ii++)

	std::fstream fs_output;
	// write function expression to file
	fs_output.open(s_FUNCTION_FILE_NAME, std::fstream::out);
	if (fs_output.is_open())
	{
		fs_output << c_tree->sGetFormula();
		fs_output.close();
	}// if (fs_output.is_open())

	 // write model data to file
	fs_output.open(s_DATA_FILE_NAME, std::fstream::out);

	bool b_success = false;
	if (fs_output.is_open())
	{
		std::string s_output_line;
		for (int ii = 0; ii < i_SAMPLE_DATA_SIZE; ii++)
		{
			for (int ij = 0; ij < i_NUMBER_OF_VARS; ij++)
			{
				fs_output << std::fixed << std::setprecision(i_WRITE_TO_FILE_PRECISION) << pd_model_data[ij][ii];
				fs_output << c_DELIMITER;
				//s_output_line += std::to_string(pd_model_data[ij][ii]);
				//s_output_line += c_DELIMITER;
			}// for (int ij = 0; ij < i_NUMBER_OF_VARS; ij++)
			 //s_output_line += std::to_string(pd_model_data[i_NUMBER_OF_VARS][ii]);
			 //fs_output << s_output_line << std::endl;
			fs_output << std::fixed << std::setprecision(i_WRITE_TO_FILE_PRECISION) << pd_model_data[i_NUMBER_OF_VARS][ii] << std::endl;
			s_output_line.clear();
		}// for (int ii = 0; ii < i_SAMPLE_DATA_SIZE; ii++)
		b_success = true;
		fs_output.close();
	}// if (fs_output.is_open())

	for (int ii = 0; ii < i_NUMBER_OF_VARS + 1; ii++)
	{
		delete[] pd_model_data[ii];
	}// for (int ii = 0; ii < i_NUMBER_OF_VARS; ii++)
	delete[] pd_model_data;
	delete c_tree;

	return b_success;
}// bool CPGAlg::bGenerateRndTreeFile(int iTreeMaxHeight)

bool CPGAlg::bGenereateFromStringTreeFile(std::string sFormula)
{
	CTree * c_tree = CTree::pcConstructFromString(sFormula);
	double ** pd_model_data = new double *[i_NUMBER_OF_VARS + 1];
	for (int ii = 0; ii < i_NUMBER_OF_VARS + 1; ii++)
	{
		pd_model_data[ii] = new double[i_SAMPLE_DATA_SIZE];
	}// for (int ii = 0; ii < i_NUMBER_OF_VARS; ii++)

	double d_rnd_value;
	std::vector<double> d_eval_data(i_NUMBER_OF_VARS);
	for (int ii = 0; ii < i_SAMPLE_DATA_SIZE; ii++)
	{
		for (int ij = 0; ij < i_NUMBER_OF_VARS; ij++)
		{
			d_rnd_value = (double)rand() / RAND_MAX;
			d_rnd_value = d_VARIABLE_MIN_VAL + d_rnd_value * (d_VARAIBLE_MAX_VAL - d_VARIABLE_MIN_VAL);
			//d_rnd_value = round(d_rnd_value * 100000) / 1000.0;  // to avoid some problembs but it doesnt really work
			pd_model_data[ij][ii] = d_rnd_value; // variable
			d_eval_data[ij] = d_rnd_value;
		}// for (int ij = 0; ij < i_NUMBER_OF_VARS; ij++)
		pd_model_data[i_NUMBER_OF_VARS][ii] = c_tree->dEvaluateFormula(d_eval_data);
	}// for (int ii = 0; ii < i_SAMPLE_DATA_SIZE; ii++)

	std::fstream fs_output;
	// write function expression to file
	fs_output.open(s_FUNCTION_FILE_NAME, std::fstream::out);
	if (fs_output.is_open())
	{
		fs_output << c_tree->sGetFormula();
		fs_output.close();
	}// if (fs_output.is_open())

	 // write model data to file
	fs_output.open(s_DATA_FILE_NAME, std::fstream::out);

	bool b_success = false;
	if (fs_output.is_open())
	{
		std::string s_output_line;
		for (int ii = 0; ii < i_SAMPLE_DATA_SIZE; ii++)
		{
			for (int ij = 0; ij < i_NUMBER_OF_VARS; ij++)
			{
				fs_output << std::fixed << std::setprecision(i_WRITE_TO_FILE_PRECISION) << pd_model_data[ij][ii];
				fs_output << c_DELIMITER;
				//s_output_line += std::to_string(pd_model_data[ij][ii]);
				//s_output_line += c_DELIMITER;
			}// for (int ij = 0; ij < i_NUMBER_OF_VARS; ij++)
			 //s_output_line += std::to_string(pd_model_data[i_NUMBER_OF_VARS][ii]);
			 //fs_output << s_output_line << std::endl;
			fs_output << std::fixed << std::setprecision(i_WRITE_TO_FILE_PRECISION) << pd_model_data[i_NUMBER_OF_VARS][ii] << std::endl;
			s_output_line.clear();
		}// for (int ii = 0; ii < i_SAMPLE_DATA_SIZE; ii++)
		b_success = true;
		fs_output.close();
	}// if (fs_output.is_open())

	for (int ii = 0; ii < i_NUMBER_OF_VARS + 1; ii++)
	{
		delete[] pd_model_data[ii];
	}// for (int ii = 0; ii < i_NUMBER_OF_VARS; ii++)
	delete[] pd_model_data;
	delete c_tree;

	return b_success;
}// bool CPGAlg::bGenereateFromStringTreeFile(std::string sFormula)

bool CPGAlg::bInitialize(CString sTest)
{
	//srand(time(NULL));
	if (!vd_data.empty())
		vd_data.clear();

	std::fstream fs_input;
	fs_input.open(sTest, std::fstream::in);
	bool b_success = false;
	if (fs_input.is_open())
	{
		std::string s_input_line;
		std::vector<double> vd_aux_vec;
		std::vector<std::string> vs_tokens;
		int i_aux_tokens_vec_size;
		while (fs_input.good())
		{
			std::getline(fs_input, s_input_line);
			if (!s_input_line.empty())
			{
				vs_tokens = vs_split_string(s_input_line, c_DELIMITER);
				i_aux_tokens_vec_size = vs_tokens.size();
				for (int ii = 0; ii < i_aux_tokens_vec_size; ii++)
				{
					vd_aux_vec.push_back(stod(vs_tokens[ii]));
				}// for (int ii = 0; ii < i_aux_tokens_vec_size; ii++)
				vd_data.push_back(vd_aux_vec);
				vd_aux_vec.clear();
			}// if (!s_input_line.empty())
		}// while (fs_input.good())
		b_success = true;
		fs_input.close();
	}// if (fs_input.is_open())
	b_initialized = b_success;

	CTree * c_tree;
	int i_initial_size = i_INITIAL_POPULATION_SIZE % 2 == 0 ? i_INITIAL_POPULATION_SIZE : i_INITIAL_POPULATION_SIZE + 1;
	int i_tree_rnd_height;
	for (int ii = 0; ii < i_initial_size; ii++)
	{
		i_tree_rnd_height = rand() % (i_TREE_MAX_HEIGHT - 1);
		c_tree = new CTree(i_tree_rnd_height);
		vc_population.push_back(c_tree);
	}// for (int ii = 0; ii < i_INITIAL_POPULATION_SIZE; ii++)
	 //CTree * tree = CTree::pcConstructFromString("+ * * 6 4 + x 7 0");
	 //vc_population[0] = tree;
	v_evaluate_trees();
	i_iteration_num = 0;
	i_iterations_without_improve = 0;
	i_modifications_iteration = 0;
	i_population_mod_iteration = 0;
	i_iteration_pop_incr = 0;
	return b_success;
}// bool CPGAlg::bInitialize(CString sTest)

void CPGAlg::vRunIteration()
{
	if (b_initialized)
	{
		i_iteration_num++;
		v_check_normalize_params();
		v_check_modifications();
		clock_t start = clock();
		v_selection();
		clock_t end = clock();
		double time_elapsed = (double)(end - start) / CLOCKS_PER_SEC;
		std::cout << "time elapsed for selection: " << time_elapsed << std::endl;
		//std::cout << "Population size: " << vc_population.size()<<std::endl;
		start = clock();
		v_delete_population();
		end = clock();
		time_elapsed = (double)(end - start) / CLOCKS_PER_SEC;
		std::cout << "time elapsed for population deletion: " << time_elapsed << std::endl;
		start = clock();
		v_parents_crossover();
		end = clock();
		time_elapsed = (double)(end - start) / CLOCKS_PER_SEC;
		std::cout << "time elapsed for parents crossover: " << time_elapsed << std::endl;
		//std::cout << "Parents size: " << vc_parents.size()<<std::endl;
		vc_parents.clear();
		start = clock();
		v_mutations();
		end = clock();
		time_elapsed = (double)(end - start) / CLOCKS_PER_SEC;
		std::cout << "time elapsed for mutation: " << time_elapsed << std::endl;
		start = clock();
		v_evaluate_trees();
		end = clock();
		time_elapsed = (double)(end - start) / CLOCKS_PER_SEC;
		std::cout << "time elapsed for evaluation: " << time_elapsed << std::endl;
		//std::cout << "Population at the end size: " << vc_population.size()<<std::endl;
	} // if (b_initialized)
}// void CPGAlg::vRunIteration()

CString CPGAlg::sGetCurrentBestTree()
{
	std::string s_formula = "";
	if (c_current_best != NULL)
		s_formula = c_current_best->sGetFormula();
	return CString(s_formula.c_str());
}// CString CPGAlg::sGetCurrentBestTree()

double CPGAlg::dGetBestDistance()
{
	return d_current_min_distance;
}// double CPGAlg::dGetBestDistance()

void CPGAlg::vPrintDataArray()
{
	if (!vd_data.empty())
	{
		int i_rows = vd_data.size();
		int i_col = vd_data[0].size();
		std::cout << "rows number: " << i_rows << ",  columns: " << i_col << std::endl << std::endl;
		for (int ii = 0; ii < i_rows; ii++)
		{
			std::cout << "row " << (ii + 1) << " :  ";
			for (int ij = 0; ij < i_col; ij++)
			{
				printf("%.8f ; ", vd_data[ii][ij]);
			}
			std::cout << std::endl;
		}
	}

}// void CPGAlg::vPrintDataArray()

std::vector<std::string> CPGAlg::vs_split_string(std::string sToSplit, char cDelimiter)
{
	std::vector<std::string> vs_tokens;
	std::string s_temp = "";
	int ii = 0;
	while (sToSplit[ii] != NULL)  // iterate throught each char
	{
		while (sToSplit[ii] != NULL && sToSplit[ii] != cDelimiter)
		{
			s_temp += sToSplit[ii];
			ii++;
		}// while (sToSplit[ii] != NULL && sToSplit[ii] != cDelimiter)
		if (s_temp[0] != NULL)
		{
			vs_tokens.push_back(s_temp);
		}// if (s_temp[0] != NULL)
		if (sToSplit[ii] == cDelimiter)
			ii++;
		s_temp.clear();
	}// while (sToSplit[ii] != NULL)

	return vs_tokens;
}// std::vector<std::string> CPGAlg::v_split_string(std::string sToSplit)

void CPGAlg::v_initialize_rnd_seed()
{
	if (!b_srand_was_called && !CTree::bGetSrandCallingStatus())
	{
		b_srand_was_called = true;
		CTree::vSetSrandCallingStatus(true);
		srand(time(NULL));
	}// if (!b_srand_was_called)
}// void CPGAlg::v_initialize_rnd_seed()

void CPGAlg::v_delete_population()
{
	int i_population_size = vc_population.size();
	for (int ii = 0; ii < i_population_size; ii++)
	{
		delete vc_population[ii];
	}// for (int ii = 0; ii < i_population_size; ii++)
}// void CPGAlg::v_delete_population()

void CPGAlg::v_delete_parents()
{
	int i_parents_size = vc_parents.size();
	for (int ii = 0; ii < i_parents_size; ii++)
	{
		delete vc_parents[ii];
	}// for (int ii = 0; ii < i_parents_size; ii++)
}// void CPGAlg::v_delete_parents()

void CPGAlg::v_delete_best_tree()
{
	if (c_current_best != NULL)
		delete c_current_best;
}// void CPGAlg::v_delete_best_tree()

void CPGAlg::v_evaluate_single_tree(CTree * c_tree)
{
	int i_data_size = vd_data.size();
	double d_distance = 0;
	double d_result;
	int i_model_nan_counter = 0;
	int i_tree_nan_counter = 0;
	int i_model_inf_counter = 0;
	int i_tree_inf_counter = 0;
	double d_parameters = vd_data[0].size();
	double d_model_result;
	double d_error;
	std::vector<double> vd_values;
	for (int ii = 0; ii < i_data_size; ii++)
	{
		vd_values = vd_data[ii];
		vd_values.resize(d_parameters - 1);
		d_result = c_tree->dEvaluateFormula(vd_values);
		d_model_result = vd_data[ii][d_parameters - 1];
		//std::cout << "tree result: " << d_result << "  ; model result: " << d_model_result << std::endl;
		//printf("tree result: %.8f  ;  model result: %.8f\n", d_result, d_model_result);
		if (!isnan(d_model_result))
		{
			if (!isnan(d_result))
			{
				if (d_model_result != INFINITY && d_model_result != -INFINITY)
				{
					if (d_result != INFINITY && d_result != -INFINITY)
					{
						d_error = (d_model_result - d_result) * (d_model_result - d_result);
						if (d_error < d_EPSILON_ACCURACY)
							d_error = 0;
						d_distance += d_error;
						//std::cout << "error function: " << d_error << std::endl;
					}// if (d_result != INFINITY && d_result != -INFINITY)
					else
					{
						i_tree_inf_counter++;
					}// else
				}// if (d_model_result != INFINITY && d_model_result != -INFINITY)
				else
				{
					i_model_inf_counter++;
					if (d_result == INFINITY || d_result == -INFINITY)
						i_tree_inf_counter++;
				}// else
			}// if (!isnan(d_result))
			else
			{
				i_tree_nan_counter++;
			}// else
		}// if (!isnan(d_model_result))
		else
		{
			i_model_nan_counter++;
			if (isnan(d_result))
				i_tree_nan_counter++;
		}// else
	}// for (int ii = 0; ii < i_data_size; ii++)
	if (i_model_nan_counter == i_tree_nan_counter && i_model_nan_counter == i_data_size)
	{
		d_distance = 0;
	}// if (i_model_nan_counter == i_tree_nan_counter && i_model_nan_counter == i_data_size)
	else if (i_model_nan_counter == i_data_size && d_distance == 0)
	{
		d_distance = INFINITY;
	}// else if (i_model_nan_counter == i_data_size && d_distance == 0)
	else if (i_tree_nan_counter == i_data_size)
	{
		d_distance = INFINITY;
	}// else if (i_tree_nan_counter == i_data_size)
	if (i_model_inf_counter == i_tree_inf_counter && i_model_inf_counter == i_data_size)
	{
		d_distance = 0;
	}// if (i_model_inf_counter == i_tree_inf_counter && i_model_inf_counter == i_data_size)
	else if (i_model_inf_counter == i_data_size && d_distance == 0)
	{
		d_distance = INFINITY;
	}// else if (i_model_inf_counter == i_data_size && d_distance == 0)
	else if (i_tree_inf_counter == i_data_size)
	{
		d_distance = INFINITY;
	}// else if (i_tree_inf_counter == i_data_size)
	c_tree->vSetDistance(d_distance);
}// void CPGAlg::v_evaluate_single_tree(CTree * c_tree)

void CPGAlg::v_evaluate_trees()
{
	int i_population_size = vc_population.size();
	double d_min_distance = INFINITY;
	int i_index_of_best = 0;
	double d_sum_of_distances = 0;
	int i_sum_of_nodes = 0;
	int i_evals_counter = 0; // debug
	for (int ii = 0; ii < i_population_size; ii++)
	{
		if (vc_population[ii]->bGetChangeStatus())  // check if tree was changed since last evaluation
		{
			v_evaluate_single_tree(vc_population[ii]);
			i_evals_counter++;
		}// if (vc_population[ii]->bGetChangeStatus())
		if (vc_population[ii]->dGetDistance() < d_min_distance)
		{
			i_index_of_best = ii;
			d_min_distance = vc_population[ii]->dGetDistance();
		}// if (vc_population[ii]->dGetDistance() < d_min_distance)
		if (!isinf(vc_population[ii]->dGetDistance()))
			d_sum_of_distances += vc_population[ii]->dGetDistance();
		i_sum_of_nodes += vc_population[ii]->iGetNumOfNodes();
	}// for (int ii = 0; ii < i_population_size; ii++)
	i_average_num_of_nodes = i_sum_of_nodes / i_population_size;
	d_average_distance = d_sum_of_distances / i_population_size;
	vi_averages_nodes_num.push_back(i_sum_of_nodes / i_population_size);
	vd_averages_distances.push_back(d_sum_of_distances / i_population_size);
	std::cout << "AVERAGE distance: " << d_average_distance << std::endl;
	std::cout << "AVERAGE num of nodes: " << i_average_num_of_nodes << std::endl;
	std::cout << "how many trees were evaluated? " << i_evals_counter << std::endl;
	std::cout << "BEST TREE OF THIS ITERATION, distance: " << d_min_distance << " ; tree: " << vc_population[i_index_of_best]->sGetFormula() << std::endl;
	if (d_min_distance < d_current_min_distance)
	{
		d_current_min_distance = d_min_distance;
		v_delete_best_tree();
		c_current_best = new CTree(*vc_population[i_index_of_best]);
		i_iterations_without_improve = 0;
	}// if (d_min_distance < d_current_min_distance)
	else
	{
		i_iterations_without_improve++;
	}// else
	vc_best_trees.push_back(new CTree(*vc_population[i_index_of_best]));
}// void CPGAlg::v_evaluate_trees()

void CPGAlg::v_selection()
{
	int i_population_size = vc_population.size();
	int i_rnd_index_1;
	int i_rnd_index_2;
	CTree * c_tree_1;
	CTree * c_tree_2;
	for (int ii = 0; ii < i_population_size; ii++)
	{
		i_rnd_index_1 = rand() % i_population_size;
		c_tree_1 = vc_population[i_rnd_index_1];
		i_rnd_index_2 = rand() % i_population_size;
		if (i_rnd_index_2 == i_rnd_index_1)
			i_rnd_index_2 = (i_rnd_index_2 + 1) % i_population_size;
		c_tree_2 = vc_population[i_rnd_index_2];
		if (c_tree_1->dGetDistance() <= c_tree_2->dGetDistance())
		{
			vc_parents.push_back(new CTree(*c_tree_1));
		}// if (c_tree_1->dGetDistance() <= c_tree_2->dGetDistance())
		else
		{
			vc_parents.push_back(new CTree(*c_tree_2));
		}// else
	}// for (int ii = 0; ii < i_population_size; ii++)
}// void CPGAlg::v_selection()

void CPGAlg::v_parents_crossover()
{
	int i_parents_size = vc_parents.size();
	int i_half_population = i_parents_size / 2;
	double d_rnd;
	CTree * c_parent_1;
	CTree * c_parent_2;
	int i_index_in_population = 0;
	for (int ii = 0; ii < i_half_population; ii++)
	{
		c_parent_1 = vc_parents[ii];
		c_parent_2 = vc_parents[i_parents_size - 1 - ii];
		d_rnd = (double)rand() / RAND_MAX;
		if (d_rnd >= 1.0 - d_cross_prob)
		{
			CTree::vCrossover(*c_parent_1, *c_parent_2);
			c_parent_1->vSetIfChanged(true);
			c_parent_2->vSetIfChanged(true);
		}// if (d_rnd >= 1.0 - d_CROSSOVER_PROBABILITY)
		else
		{
			c_parent_1->vSetIfChanged(false);
			c_parent_2->vSetIfChanged(false);
		}// else
		vc_population[i_index_in_population] = c_parent_1;
		i_index_in_population++;
		vc_population[i_index_in_population] = c_parent_2;
		i_index_in_population++;
	}// for (int ii = 0; ii < i_half_population; ii++)
}// void CPGAlg::v_parents_crossover()

void CPGAlg::v_mutations()
{
	int i_population_size = vc_population.size();
	double d_rnd;
	for (int ii = 0; ii < i_population_size; ii++)
	{
		d_rnd = (double)rand() / RAND_MAX;
		if (d_rnd >= 1.0 - d_mut_prob)
		{
			vc_population[ii]->vMutation();
			vc_population[ii]->vSetIfChanged(true);
		}// if (d_rnd >= 1.0 - d_MUTATION_PROBABILITY)
		else
		{
			vc_population[ii]->vSetIfChanged((vc_population[ii]->bGetChangeStatus() || false));  // leave true status (changed) if tree was created as crossover effect
		}// else
	}// for (int ii = 0; ii < i_population_size; ii++)
}// void CPGAlg::v_mutations()

void CPGAlg::v_mixin_bests()
{
	int i_population_size = vc_population.size();
	std::vector<std::pair<double, int>> v_population_distance_index;
	for (int ii = 0; ii < i_population_size; ii++)
	{
		v_population_distance_index.push_back(std::pair<double, int>(vc_population[ii]->dGetDistance(), ii));
	}// for (int ii = 0; ii < i_population_size; ii++)
	std::sort(v_population_distance_index.begin(), v_population_distance_index.end(), std::greater<std::pair<double, int>>());
	std::vector<std::pair<double, int>> v_bests_distance_index;
	int i_num_of_bests = vc_best_trees.size();
	for (int ii = 0; ii < i_num_of_bests; ii++)
	{
		v_bests_distance_index.push_back(std::pair<double, int>(vc_best_trees[ii]->dGetDistance(), ii));
	}// for (int ii = 0; ii < i_population_size; ii++)
	std::sort(v_bests_distance_index.begin(), v_bests_distance_index.end());
	int i_index_to_replace;
	int i_index_replacing;
	for (int ii = 0; ii < i_HOW_MANY_TO_MIXIN && ii < i_num_of_bests; ii++)
	{
		i_index_to_replace = v_population_distance_index[ii].second;
		delete vc_population[i_index_to_replace];
		i_index_replacing = v_bests_distance_index[ii].second;
		vc_population[i_index_to_replace] = new CTree(*vc_best_trees[i_index_replacing]);
	}// for (int ii = 0; ii < i_HOW_MANY_TO_MIXIN && ii < i_num_of_bests; ii++)
}// void CPGAlg::v_mixin_bests()

void CPGAlg::v_check_modifications()
{
	if (i_iteration_num % i_MIX_IN_BESTS_INTERVAL == 0)
	{
		v_mixin_bests();
		v_mixin_mutated_best_tree_family();
	}// if (i_iteration_num == i_MIX_IN_BESTS_INTERVAL)
	if (i_iterations_without_improve >= i_ACCEPTABLE_ITERATIONS_WITHOUT_IMPROVEMENT)
	{
		v_increase_probabilities();
		i_modifications_iteration = i_iteration_num;
	}// if (i_iterations_without_improve >= i_ACCEPTABLE_ITERATIONS_WITHOUT_IMPROVEMENT)
	if (i_iteration_num % i_LINEAR_REGRESSION_SAMPLE == 0 && b_nodes_num_optimization_needed())
	{
		v_replace_worsts_with_new();
		i_population_mod_iteration = i_iteration_num;
	}// if (i_iteration_num % i_LINEAR_REGRESSION_SAMPLE == 0 && b_nodes_num_optimization_needed())
	if (i_average_num_of_nodes > i_ACCEPTABLE_NUM_OF_NODES && i_iteration_num - i_population_mod_iteration >= i_POPULATION_MOD_INTERVAL)
	{
		v_replace_worsts_with_new();
		i_population_mod_iteration = i_iteration_num;
	}// if (i_average_num_of_nodes > i_ACCEPTABLE_NUM_OF_NODES && i_iteration_num - i_population_mod_iteration >= i_POPULATION_MOD_INTERVAL)
	if (i_iterations_without_improve >= i_ACCEPTABLE_ITERATIONS_WITHOUT_IMPROVEMENT && i_iteration_num - i_population_mod_iteration >= i_POPULATION_MOD_INTERVAL)
	{
		v_replace_worsts_with_new();
		i_population_mod_iteration = i_iteration_num;
	}// if (i_iterations_without_improve >= i_ACCEPTABLE_ITERATIONS_WITHOUT_IMPROVEMENT && i_iteration_num - i_population_mod_iteration >= i_POPULATION_MOD_INTERVAL)
	if (i_average_num_of_nodes < i_ACCEPTABLE_NUM_OF_NODES * 2 && i_iterations_without_improve >= i_ACCEPTABLE_ITERATIONS_WITHOUT_IMPROVEMENT && !b_population_increased)
	{
		b_population_increased = true;
		i_iteration_pop_incr = i_iteration_num;
		v_increase_population_size();
	}// if (i_iterations_without_improve >= i_ACCEPTABLE_ITERATIONS_WITHOUT_IMPROVEMENT && !b_population_increased)
}// void CPGAlg::v_check_modifications()

void CPGAlg::v_check_normalize_params()
{
	if (i_iteration_num > i_modifications_iteration)
	{
		v_normalize_probabilities();
	}// if (i_iteration_num > i_modifications_iteration)
	if (i_iteration_num - i_iteration_pop_incr > i_ITERATIONS_WITH_INCREASED_POPULATION && b_population_increased)
	{
		v_normalize_population_size();
		b_population_increased = false;
	}// if (i_iterations_without_improve >= i_ACCEPTABLE_ITERATIONS_WITHOUT_IMPROVEMENT && i_iteration_num - i_iteration_pop_incr > i_ITERATIONS_WITH_INCREASED_POPULATION && b_population_increased)
}// void CPGAlg::v_check_normalize_params()

void CPGAlg::v_increase_probabilities()
{
	d_mut_prob = d_INCREASED_MUTATION_PROB;
	d_cross_prob = d_INCREASED_CROSSOVER_PROB;
}// void CPGAlg::v_increase_probabilities()

void CPGAlg::v_normalize_probabilities()
{
	d_mut_prob = d_MUTATION_PROBABILITY;
	d_cross_prob = d_CROSSOVER_PROBABILITY;
}// void CPGAlg::v_normalize_probabilities()

bool CPGAlg::b_nodes_num_optimization_needed()
{
	// y = ax + b, where x: average number of nodes at i-th iteration, y: average distance at i-th iteration
	// a parameter, the only one we need cause its tangent of line's angle of slope
	double d_a;
	double d_xy_sum = 0;
	double d_x_sum = 0;
	double d_y_sum = 0;
	double d_x_sqr_sum = 0;
	for (int ii = 0; ii < i_LINEAR_REGRESSION_SAMPLE; ii++)
	{
		d_xy_sum += vi_averages_nodes_num[ii] * vd_averages_distances[ii];
		d_x_sum += vi_averages_nodes_num[ii];
		d_y_sum += vd_averages_distances[ii];
		d_x_sqr_sum += vi_averages_nodes_num[ii] * vi_averages_nodes_num[ii];
	}// for (int ii = i_start_from; ii < i_iteration_num; ii++)
	double d_sqr_x_sum = d_x_sum * d_x_sum;
	d_a = (i_LINEAR_REGRESSION_SAMPLE * d_xy_sum - d_x_sum * d_y_sum) / (i_LINEAR_REGRESSION_SAMPLE * d_x_sqr_sum - d_sqr_x_sum);

	bool b_optimize_num_of_nodes = false;
	if (d_a > 0)  // increasing function
	{
		b_optimize_num_of_nodes = true;
	}// if (d_a > 0)
	vi_averages_nodes_num.clear();
	vd_averages_distances.clear();

	return b_optimize_num_of_nodes;
}// // bool CPGAlg::b_nodes_num_optimization_needed()

void CPGAlg::v_replace_worsts_with_new()
{
	int i_population_size = vc_population.size();
	std::vector<std::pair<double, int>> v_population_distance_index;
	for (int ii = 0; ii < i_population_size; ii++)
	{
		v_population_distance_index.push_back(std::pair<double, int>(vc_population[ii]->dGetDistance(), ii));
	}// for (int ii = 0; ii < i_population_size; ii++)
	std::sort(v_population_distance_index.begin(), v_population_distance_index.end(), std::greater<std::pair<double, int>>());
	int i_index_to_replace;
	int i_rnd_height;
	for (int ii = 0; ii < i_HOW_MANY_TO_REPLACE_WTH_NEW; ii++)
	{
		i_index_to_replace = v_population_distance_index[ii].second;
		delete vc_population[i_index_to_replace];

		if (i_average_num_of_nodes <= i_ACCEPTABLE_NUM_OF_NODES)
			i_rnd_height = rand() % (i_TREE_MAX_HEIGHT + 1);
		else if (i_average_num_of_nodes > i_ACCEPTABLE_NUM_OF_NODES && i_average_num_of_nodes <= i_ACCEPTABLE_NUM_OF_NODES * 3)
			i_rnd_height = rand() % (i_TREE_MAX_HEIGHT);
		else
			i_rnd_height = rand() % (i_TREE_MAX_HEIGHT - 1);

		vc_population[i_index_to_replace] = new CTree(i_rnd_height);
		v_evaluate_single_tree(vc_population[i_index_to_replace]);
	}// for (int ii = 0; ii < i_HOW_MANY_TO_REPLACE_WTH_NEW; ii++)
}// void CPGAlg::v_replace_worsts_with_new()

void CPGAlg::v_increase_population_size()
{
	int i_population_size = vc_population.size();
	int i_new_size = i_population_size * d_INCREASE_POPULATION_FACTOR;
	if (i_new_size % 2 != 0)
		i_new_size++;
	int i_difference = i_new_size - i_population_size;
	int i_rnd_height;
	CTree * c_tree;
	for (int ii = 0; ii < i_difference; ii++)
	{
		if (i_average_num_of_nodes > i_ACCEPTABLE_NUM_OF_NODES)
			i_rnd_height = rand() % (i_TREE_MAX_HEIGHT - 1);
		else
			i_rnd_height = rand() % (i_TREE_MAX_HEIGHT);
		c_tree = new CTree(i_rnd_height);
		v_evaluate_single_tree(c_tree);
		vc_population.push_back(c_tree);
	}// for (int ii = 0; ii < i_difference; ii++)
}// void CPGAlg::v_increase_population_size()

void CPGAlg::v_normalize_population_size()
{
	int i_current_population_size = vc_population.size();
	std::vector<std::pair<double, int>> v_population_distance_index;
	for (int ii = 0; ii < i_current_population_size; ii++)
	{
		v_population_distance_index.push_back(std::pair<double, int>(vc_population[ii]->dGetDistance(), ii));
	}// for (int ii = 0; ii < i_population_size; ii++)
	std::sort(v_population_distance_index.begin(), v_population_distance_index.end(), std::greater<std::pair<double, int>>());
	int i_how_many_to_delete = i_current_population_size - i_INITIAL_POPULATION_SIZE;
	int i_index_to_delete;
	std::vector<CTree *> vc_aux_population;
	std::vector<int> vi_to_delete;
	for (int ii = 0; ii < i_how_many_to_delete; ii++)
	{
		i_index_to_delete = v_population_distance_index[ii].second;
		vi_to_delete.push_back(i_index_to_delete);
	}// for (int ii = 0; ii < i_how_many_to_delete; ii++)
	std::sort(vi_to_delete.begin(), vi_to_delete.end());
	int ij = 0;
	for (int ii = 0; ii < i_current_population_size; ii++)
	{
		if (ij >= i_how_many_to_delete || ii != vi_to_delete[ij])
		{
			vc_aux_population.push_back(vc_population[ii]);
		}// if (ii != vi_to_delete[ij])
		else
		{
			ij++;
		}// else
	}// for (int ii = 0; ii < i_how_many_to_delete; ii++)
	for (int ii = 0; ii < i_how_many_to_delete; ii++)
	{
		delete vc_population[vi_to_delete[ii]];
	}// for (int ii = 0; ii < i_how_many_to_delete; ii++)
	vc_population = vc_aux_population;
}// void CPGAlg::v_normalize_population_size()

void CPGAlg::v_mixin_mutated_best_tree_family()
{
	if (c_current_best != NULL)
	{
		std::vector<CTree*> vc_best_tree_family;
		for (int ii = 0; ii < i_BEST_TREE_FAMILY_SIZE; ii++)
		{
			vc_best_tree_family.push_back(new CTree(*c_current_best));
			vc_best_tree_family[ii]->vMutation();
			v_evaluate_single_tree(vc_best_tree_family[ii]);
		}// for (int ii = 0; ii < i_BEST_TREE_FAMILY_SIZE; ii++)

		int i_population_size = vc_population.size();
		std::vector<std::pair<double, int>> v_population_distance_index;
		for (int ii = 0; ii < i_population_size; ii++)
		{
			v_population_distance_index.push_back(std::pair<double, int>(vc_population[ii]->dGetDistance(), ii));
		}// for (int ii = 0; ii < i_population_size; ii++)
		std::sort(v_population_distance_index.begin(), v_population_distance_index.end(), std::greater<std::pair<double, int>>());
		std::vector<std::pair<double, int>> v_bests_distance_index;
		for (int ii = 0; ii < i_BEST_TREE_FAMILY_SIZE; ii++)
		{
			v_bests_distance_index.push_back(std::pair<double, int>(vc_best_tree_family[ii]->dGetDistance(), ii));
		}// for (int ii = 0; ii < i_population_size; ii++)
		std::sort(v_bests_distance_index.begin(), v_bests_distance_index.end());
		int i_index_to_replace;
		int i_index_replacing;
		for (int ii = 0; ii < i_BEST_TREE_FAMILY_SIZE; ii++)
		{
			i_index_to_replace = v_population_distance_index[ii].second;
			delete vc_population[i_index_to_replace];
			i_index_replacing = v_bests_distance_index[ii].second;
			vc_population[i_index_to_replace] = new CTree(*vc_best_tree_family[i_index_replacing]);
		}// for (int ii = 0; ii < i_HOW_MANY_TO_MIXIN && ii < i_num_of_bests; ii++)
		for (int ii = 0; ii < i_BEST_TREE_FAMILY_SIZE; ii++)
		{
			delete vc_best_tree_family[ii];
		}// for (int ii = 0; ii < i_BEST_TREE_FAMILY_SIZE; ii++)
		vc_best_tree_family.clear();
	}// if (c_current_best != NULL)
}// void CPGAlg::v_mixin_mutated_best_tree_family()

bool CPGAlg::b_srand_was_called = false;
