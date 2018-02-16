#include "CTree.h"
#include <iostream>
// CTree part

CTree::CTree()
{
	pc_root = nullptr;
	b_changed_since_eval = true;
	d_distance_measure = INFINITY;
}// CTree::CTree()

CTree::CTree(std::string sInputFormula)
{
	b_changed_since_eval = true;
	d_distance_measure = INFINITY;
	int *i_current_pos = new int(0);
	pc_root = new CNode(nullptr);
	i_num_of_nodes = 0;
	v_construct_node(pc_root, sInputFormula, i_current_pos);
	delete i_current_pos;
}// CTree::CTree(std::string sInputFomrula)

CTree::CTree(CTree const & cOtherTree)
{
	pc_root = new CNode(nullptr);
	i_num_of_nodes = cOtherTree.i_num_of_nodes;
	d_distance_measure = cOtherTree.d_distance_measure;
	b_changed_since_eval = cOtherTree.b_changed_since_eval;
	v_copy_node(pc_root, *cOtherTree.pc_root);
}// CTree::CTree(CTree const & cOtherTree)

CTree::CTree(int iMaxHeight)
{
	if (iMaxHeight < 0)  // max height can't be negative, don't throw exception or pass error code, just gemerate tree with height = 0 - just non-operator root
		iMaxHeight = 0;
	v_initialize_rnd_seed();
	pc_root = pc_gen_rnd_root(iMaxHeight);
	i_num_of_nodes = 1;
	b_changed_since_eval = true;
	d_distance_measure = INFINITY;
	int i_root_children_num;
	if (b_token_is_operator(pc_root->s_value))
	{
		i_root_children_num = msi_OPERATOR_CHILDREN_NUM.at(pc_root->s_value);
	}// if (b_token_is_operator(pc_root->s_value))
	else
	{
		i_root_children_num = 0;
		vc_leaves.push_back(pc_root);
		if (pc_root->b_is_var())
			msd_vars.insert(std::pair<std::string, double>(pc_root->s_value, 0));
	}// else
	v_gen_rnd_tree(iMaxHeight, i_root_children_num);
	//v_make_vars_and_leaves();
}// CTree::CTree(int iMaxHeight)

CTree::~CTree()
{
	msd_vars.clear();
	vc_leaves.clear();
	delete pc_root;
}// CTree::~CTree()

CTree * CTree::pcConstructFromString(std::string sFormula)
{
	std::string s_fixed = s_fix_input_formula(sFormula);
	if (s_fixed == sUNALLOWED_EXP)
		return nullptr;

	return new CTree(s_fixed);
}// CTree * CTree::pcConstructFromString(std::string sFormula)

void CTree::v_construct_node(CNode * pcNode, std::string const & sFormula, int *piCurrentPos)
{
	i_num_of_nodes++;
	std::string s_token = sNextToken(sFormula, piCurrentPos);
	pcNode->s_value = s_token;
	int i_children_num;
	if (b_token_is_operator(pcNode->s_value))
	{
		i_children_num = msi_OPERATOR_CHILDREN_NUM.at(pcNode->s_value);
	}// if (b_token_is_operator(pcNode->s_value))
	else
	{
		i_children_num = 0;
		vc_leaves.push_back(pcNode);
		if (pcNode->b_is_var())
			msd_vars.insert(std::pair<std::string, double>(pcNode->s_value, 0));
	}// else
	CNode* pc_child;
	for (int ii = 0; ii < i_children_num; ii++)
	{
		pc_child = new CNode(pcNode);
		v_construct_node(pc_child, sFormula, piCurrentPos);
		pcNode->vc_children.push_back(pc_child);
	}// for (int ii = 0; ii < i_children_num; ii++)
}// void CTree::v_construct_node(CNode * pcNode, std::string const & sFormula, int *piCurrentPos)

void CTree::v_copy_node(CNode * pcNode, CNode & cOtherNode)
{
	pcNode->s_value = cOtherNode.s_value;
	int i_children_num = cOtherNode.vc_children.size();
	if (i_children_num == 0)  // leaf case
	{
		vc_leaves.push_back(pcNode);
		if (pcNode->b_is_var())
			msd_vars.insert(std::pair<std::string, double>(pcNode->s_value, 0));
	}// if (i_children_num == 0)
	CNode * pc_child;
	for (int ii = 0; ii < i_children_num; ii++)
	{
		pc_child = new CNode(pcNode);
		v_copy_node(pc_child, *cOtherNode.vc_children[ii]);
		pcNode->vc_children.push_back(pc_child);
	}// for (int ii = 0; ii < i_children_num; ii++)
}// void CTree::v_copy_node(CNode * pcNode, CNode & pcOtherNode)

std::string CTree::sGetFormula()
{
	return s_aux_preorder_walk(pc_root);
}// std::string CTree::sGetFormula()

std::string CTree::sGetVars()
{
	std::string s_result = sEMPTY;
	for (std::map <std::string, double > ::iterator it = msd_vars.begin(); it != msd_vars.end(); it++)
	{
		s_result += it->first+cSPACE;
	}// for (std::map <std::string, double > ::iterator it = msd_vars.begin(); it != msd_vars.end(); it++)
	return s_result;
}// std::string CTree::sGetVars()

std::string CTree::sGetLeaves()
{
	std::string s_result = sEMPTY;
	int i_leaves_num = vc_leaves.size();
	for (int ii = 0; ii < i_leaves_num; ii++)
	{
		s_result += vc_leaves[ii]->s_value+" ";
	}// for (int ii = 0; ii < i_vars_num; ii++)
	return s_result;
}

double CTree::dEvaluateFormula(std::vector<double> const &vdValues)
{
	int i_values_num = vdValues.size();
	int i_var_num = msd_vars.size();
	double d_result = -1;
	if (i_values_num < i_var_num)
	{
		err_code = CTreeError::ERR_NOT_ENOUGH_VALUES;
	}// else if (i_values_num < i_var_num)
	else
	{
		int i_index = 0;
		for (std::map<std::string, double>::iterator it = msd_vars.begin(); it != msd_vars.end(); it++)
		{
			it->second = vdValues[i_index];
			i_index++;
		}// for (std::map<std::string, double>::iterator it = msd_vars.begin(); it != msd_vars.end(); it++)
		d_result= d_aux_evaluate(pc_root);
		err_code = CTreeError::ERR_OK;
	}// else
	return d_result;
}// double CTree::dEvaluateFormula(std::vector<double> const &vdValues)

CTree CTree::operator+(CTree const & cOtherTree)
{
	//srand(time(NULL));
	// creating copy of this tree
	CTree c_new_tree(*this);
	int i_leaves_num = c_new_tree.vc_leaves.size();
	int i_which_leaf = rand() % i_leaves_num;
	CNode *&pc_selected = c_new_tree.vc_leaves[i_which_leaf];
	// we have to find which child is pc_selected for its parent
	int i_parent_children_num = pc_selected->pc_parent->vc_children.size();
	int i_child_num = 0;
	bool bf_found_child = false;
	for (int ii = 0; ii < i_parent_children_num && !bf_found_child; ii++)
	{
		if (pc_selected->pc_parent->vc_children[ii] == pc_selected)
		{
			i_child_num = ii;
			bf_found_child = true;
		}// if (pc_selected->pc_parent->vc_children[ii] == pc_selected)
	}// for (int ii = 0; ii < i_parent_children_num && !bf_found_child; ii++)
	CNode * pc_parent_of_selected = pc_selected->pc_parent;
	delete pc_selected;
	pc_selected = new CNode(pc_parent_of_selected);
	v_copy_node(pc_selected, *cOtherTree.pc_root);
	pc_selected->pc_parent->vc_children[i_child_num] = pc_selected;
	return c_new_tree;
}// CTree CTree::operator+(CTree const & cOtherTree)

CTree & CTree::operator=(CTree const & cOtherTree)
{
	msd_vars.clear();
	vc_leaves.clear();
	CNode *pc_new_root = new CNode(nullptr);
	d_distance_measure = cOtherTree.d_distance_measure;
	b_changed_since_eval = cOtherTree.b_changed_since_eval;
	i_num_of_nodes = cOtherTree.i_num_of_nodes;
	v_copy_node(pc_new_root, *cOtherTree.pc_root);
	if (pc_root != nullptr)
		delete pc_root;
	pc_root = pc_new_root;
	return (*this);
}// CTree & CTree::operator=(CTree const & cOtherTree)

CTreeError CTree::errGetErrorCode()
{
	return err_code;
}// CTreeError CTree::errGetErrorCode()

void CTree::vMutation()
{
	CNode** pc_rnd = pc_get_rnd_node(true);
	double d_rnd = c_UNIF_REAL_PROB_DIST(c_GENERATOR);
	
	if (d_rnd >= 1.0 - d_SELF_NODE_SWAP_PROB && *pc_rnd != pc_root)  // replace rnd node with another self rnd node
	{
		v_swap_with_self_node(*pc_rnd);
	}// if (d_rnd >= 1.0 - d_SELF_NODE_SWAP_PROB && *pc_rnd != pc_root)
	else  // replace random node with new rnd tree (node)
	{
		v_replace_with_new_node(*pc_rnd);
	}// else
	
	//std::cout << "selected: " << (*pc_rnd)->s_value << std::endl;  //debug
	v_make_vars_and_leaves();
}// void CTree::vMutation()

void CTree::vCrossover(CTree & cTree1, CTree & cTree2)
{
	if (cTree1.i_num_of_nodes > 1 && cTree2.i_num_of_nodes > 1)
	{
		CNode * pc_split_node_tree1 = *cTree1.pc_get_rnd_node(false);
		CNode * pc_split_node_tree2 = *cTree2.pc_get_rnd_node(false);
		CNode * pc_parent1 = pc_split_node_tree1->pc_parent;
		CNode * pc_parent2 = pc_split_node_tree2->pc_parent;
		int i_child_num1 = pc_parent1->vc_children.size();
		int i_which_child1 = -1;
		for (int ii = 0; ii < i_child_num1 && i_which_child1 == -1; ii++)
		{
			if (pc_parent1->vc_children[ii] == pc_split_node_tree1)
			{
				i_which_child1 = ii;
			}// if (pc_parent1->vc_children[ii] == pc_split_node_tree1)
		}// for (int ii = 0; ii < i_child_num1 && i_which_child1 == -1; ii++)

		int i_child_num2 = pc_parent2->vc_children.size();
		int i_which_child2 = -1;
		for (int ii = 0; ii < i_child_num2 && i_which_child2 == -1; ii++)
		{
			if (pc_parent2->vc_children[ii] == pc_split_node_tree2)
			{
				i_which_child2 = ii;
			}// if (pc_parent2->vc_children[ii] == pc_split_node_tree2)
		}// for (int ii = 0; ii < i_child_num2 && i_which_child2 == -1; ii++)

		pc_parent1->vc_children[i_which_child1] = pc_split_node_tree2;
		pc_split_node_tree2->pc_parent = pc_parent1;

		pc_parent2->vc_children[i_which_child2] = pc_split_node_tree1;
		pc_split_node_tree1->pc_parent = pc_parent2;

		cTree1.v_make_vars_and_leaves();
		cTree2.v_make_vars_and_leaves();
	}// if (cTree1.i_num_of_nodes > 1 && cTree2.i_num_of_nodes > 1)
	
}// void CTree::vCrossover(CTree & cTree1, CTree & cTree2)

double CTree::dGetDistance()
{
	return d_distance_measure;
}// double CTree::dGetDistance()

void CTree::vSetDistance(double dDistance)
{
	d_distance_measure = dDistance;
}// void CTree::vSetDistance()

bool CTree::bGetChangeStatus()
{
	return b_changed_since_eval;
}// bool CTree::bGetChangeStatus()

void CTree::vSetIfChanged(bool b_changed)
{
	b_changed_since_eval = b_changed;
}// void CTree::vSetIfChanged(bool b_changed)

int CTree::iGetNumOfNodes()
{
	return i_num_of_nodes;
}// int CTree::iGetNumOfNodes()

bool CTree::bGetSrandCallingStatus()
{
	return b_srand_was_called;
}// bool CTree::bGetSrandCallingStatus()

void CTree::vSetSrandCallingStatus(bool bStatus)
{
	b_srand_was_called = bStatus;
}// void CTree::vSetSrandCallingStatus(bool bStatus)

std::string CTree::s_fix_input_formula(std::string sFormula)
{
	std::string s_fixed = sEMPTY;
	std::string s_token;
	int i_current_node_need = 0;
	bool bf_fix_need = true;
	int i_current_pos = 0;
	int *pi_current_pos = &i_current_pos;
	s_token = sNextToken(sFormula, pi_current_pos);
	i_current_pos = 0;
	if (!b_token_is_operator(s_token))
	{
		sFormula = cSPACE + sFormula;
		sFormula = ADDITION + sFormula;
	}// if (!b_token_is_operator(sToken))
	int i_formula_len = sFormula.length();
	int i_num_of_tokens = iNumOfTokensInString(sFormula);
	for (int ii = 0; ii < i_num_of_tokens && bf_fix_need; ii++)
	{
		s_token = sNextToken(sFormula, pi_current_pos);
	
		if (s_token == ADDITION || s_token == SUBTRACTION || s_token == MULTIPLICATION || s_token == DIVISION)
		{
			if (i_current_node_need > 0)
			{
				i_current_node_need--;
			}// if (i_current_node_need > 0)
			i_current_node_need += 2;
			s_fixed += s_token;
		}// if (sToken == ADDITION || sToken = SUBTRACTION || sToken = MULTIPLICATION || sToken == DIVISION)
		else if (s_token == SIN || s_token == COS)
		{
			if (i_current_node_need > 0)
			{
				i_current_node_need--;
			}// if (i_current_node_need > 0)
			i_current_node_need += 1;
			s_fixed += s_token;
		}// else if (sToken == SIN || sToken == COS)
		else
		{
			if (i_current_node_need == 0)
			{
				bf_fix_need = false;
			}// if (i_current_node_need == 0)
			else
			{
				if (b_token_is_var(s_token))
				{
					i_current_node_need--;
					s_fixed += s_token;
				}// if (b_token_is_var(sToken))
				else if (bIsInteger(s_token))
				{
					i_current_node_need--;
					s_fixed += s_token;
				}// else if (bIsInteger(sToken))
				else
				{
					i_current_node_need = 0;
					s_fixed = sUNALLOWED_EXP;
					bf_fix_need = false;
				}// else
			}// else
			if (i_current_node_need == 0)
			{
				bf_fix_need = false;
			}// if (i_current_node_need == 0)
		}// else
		s_fixed += cSPACE;
	}// for (int ii = 0; ii < i_formula_len && bf_fix_need; ii++)
	s_fixed=s_fixed.substr(0, s_fixed.length() - 1);  // to trim ending space
	for (int ii = 0; ii < i_current_node_need; ii++)
	{
		s_fixed += cSPACE;
		s_fixed += cFIX_NUM;
	}// for (int ii = 0; ii < i_current_node_need; ii++)
	return s_fixed;
}// std::string CTree::s_fix_input_formula(std::string sFormula)

std::string CTree::s_aux_preorder_walk(CNode * pcNode)
{
	std::string s_result = pcNode->s_value + cSPACE;
	int i_children_num = pcNode->vc_children.size();
	for (int ii = 0; ii < i_children_num; ii++)
	{
		s_result += s_aux_preorder_walk(pcNode->vc_children[ii]);
	}// for (int ii = 0; ii < i_children_num; ii++)
	return s_result;
}// std::string CTree::s_aux_preorder_walk(CNode * c_node)

double CTree::d_aux_evaluate(CNode * pcNode)
{
	if (pcNode->s_value == ADDITION)
	{
		return d_aux_evaluate(pcNode->vc_children[0]) + d_aux_evaluate(pcNode->vc_children[1]);
	}// if (pcNode->s_value == ADDITION)
	else if (pcNode->s_value == SUBTRACTION)
	{
		return d_aux_evaluate(pcNode->vc_children[0]) - d_aux_evaluate(pcNode->vc_children[1]);
	}// else if (pcNode->s_value == SUBTRACTION)
	else if (pcNode->s_value == MULTIPLICATION)
	{
		return d_aux_evaluate(pcNode->vc_children[0]) * d_aux_evaluate(pcNode->vc_children[1]);
	}// else if (pcNode->s_value == MULTIPLICATION)
	else if (pcNode->s_value == DIVISION)
	{
		return d_aux_evaluate(pcNode->vc_children[0]) / d_aux_evaluate(pcNode->vc_children[1]);
	}// else if (pcNode->s_value == DIVISION)
	else if (pcNode->s_value == SIN)
	{
		return sin(d_aux_evaluate(pcNode->vc_children[0]));
	}// else if (pcNode->s_value == SIN)
	else if (pcNode->s_value == COS)
	{
		return cos(d_aux_evaluate(pcNode->vc_children[0]));
	}// else if (pcNode->s_value == COS)
	else
	{
		if (pcNode->b_is_var())  // check if node is variable, cause then we have to use value from vc_vars vector
		{
			return d_get_value_of_var(pcNode->s_value);
		}// if (pcNode->b_is_var())
		else
		{
			return stod(pcNode->s_value);
		}// else
	}// else
}// double CTree::d_aux_evaluate(CNode * pcNode)

double CTree::d_get_value_of_var(std::string const & sVarName)
{
	return msd_vars.at(sVarName);
}// double CTree::d_get_value_of_var(std::string const & sVarName)

CNode * CTree::pc_gen_rnd_root(int iMaxHeight)
{
	double d_rnd = c_UNIF_REAL_PROB_DIST(c_GENERATOR);
	CNode *pc_generated_root = new CNode();
	if (iMaxHeight > 0 && d_rnd >= 1.0 - d_ROOT_OPERATOR_PROBABILITY)  // 90% chance to generate operator
	{
		pc_generated_root->s_value = s_gen_rnd_operator();
	}// if (iMaxheight > 0 && d_rnd >= 1.0 - 0.9)
	else
	{
		pc_generated_root->s_value = s_gen_rnd_var_or_val();
	}// else
	return pc_generated_root;
}// CNode * CTree::pc_gen_rnd_root(int iMaxHeight)

void CTree::v_gen_rnd_tree(int iMaxHeight, int iRootChildrenNum)
{
	CNode* pc_child;
	for (int ii = 0; ii < iRootChildrenNum; ii++)
	{
		pc_child = new CNode(iMaxHeight, pc_root);
		v_gen_rnd_node(iMaxHeight, pc_child);
		pc_root->vc_children.push_back(pc_child);
	}// for (int ii = 0; ii < iRootChildrenNum; ii++)
}// void CTree::v_gen_rnd_tree(int iMaxHeight)

void CTree::v_gen_rnd_node(int iMaxHeight, CNode * pcNode)
{
	i_num_of_nodes++;
	if (pcNode->i_depth < iMaxHeight)  // operator node can be generated
	{
		double d_rnd = c_UNIF_REAL_PROB_DIST(c_GENERATOR);
		if (d_rnd >= 1.0 - d_OPERATOR_GENERATE_PROBABILITY)  // 60% chance to generate operator
		{
			pcNode->s_value = s_gen_rnd_operator();
			int i_children_num = msi_OPERATOR_CHILDREN_NUM.at(pcNode->s_value);
			CNode* pc_child;
			for (int ii = 0; ii < i_children_num; ii++)
			{
				pc_child = new CNode(iMaxHeight, pcNode);
				v_gen_rnd_node(iMaxHeight, pc_child);
				pcNode->vc_children.push_back(pc_child);
			}// for (int ii = 0; ii < i_children_num; ii++)
		}// if (d_rnd >= 1.0 - d_OPERATOR_GENERATE_PROBABILITY)
		else
		{
			pcNode->s_value = s_gen_rnd_var_or_val();
			vc_leaves.push_back(pcNode);
			if(pcNode->b_is_var())
				msd_vars.insert(std::pair<std::string, double>(pcNode->s_value, 0));
		}//else		
	}
	else  // maxHeight has been reached, only var or val can be generated
	{
		pcNode->s_value = s_gen_rnd_var_or_val();
		vc_leaves.push_back(pcNode);
		if (pcNode->b_is_var())
			msd_vars.insert(std::pair<std::string, double>(pcNode->s_value, 0));
	}// else
}// void CTree::v_gen_rnd_node(int iMaxHeight)

CNode ** CTree::pc_get_rnd_node(bool bfRootIsAcceptable)
{
	int i_rnd;
	if (bfRootIsAcceptable)
		i_rnd = rand() % i_num_of_nodes;
	else
		i_rnd = rand() % (i_num_of_nodes - 1) + 1;
	//std::cout << "\nrandomly selected: " << i_rnd << std::endl;  //debug
	return pc_breadth_select(i_rnd);
}// CNode * CTree::pc_get_rnd_node()

void CTree::v_replace_with_new_node(CNode *& pcSelected)
{
	int i_new_tree_height;
	if (pcSelected->i_depth > i_MAX_HEIGHT)
		i_new_tree_height = 0;
	else
		i_new_tree_height = rand() % i_MAX_HEIGHT;  // so upper bound is max height -1 !!!
	CTree * pc_new_tree = new CTree(i_new_tree_height);  // create rnd tree with rnd height
	//std::cout << "new_tree: " << pc_new_tree->sGetFormula() << std::endl;  //debug
	CNode * pc_selected_parent = pcSelected->pc_parent;
	delete pcSelected;
	pcSelected = pc_new_tree->pc_root;
	pcSelected->pc_parent = pc_selected_parent;
}// void CTree::v_replace_with_new_node(CNode * pcSelected)

void CTree::v_swap_with_self_node(CNode *& pcSelected)
{
	bool bf_go_up = true;
	CNode ** pc_aux_walker = &pcSelected;
	CNode * pc_came_from_child;
	double d_rnd;
	do
	{
		pc_came_from_child = *pc_aux_walker;
		pc_aux_walker = &(*pc_aux_walker)->pc_parent;
		d_rnd = c_UNIF_REAL_PROB_DIST(c_GENERATOR);
		if (d_rnd < 1.0 - d_GO_UP_DOWN_PROB)
		{
			bf_go_up = false;
		}// if (d_rnd < 1.0 - d_GO_UP_DOWN_PROB)
		if ((*pc_aux_walker)->vc_children.size() == 1)
		{
			bf_go_up = true;
		}// if (pc_aux_walker->vc_children.size() == 1)
		if ((*pc_aux_walker)->pc_parent == nullptr)
		{
			bf_go_up = false;
		}// if (pc_aux_walker->pc_parent == nullptr)

	} while (bf_go_up);
	// do while (bf_go_up)

	int i_walker_children_num = (*pc_aux_walker)->vc_children.size();
	if (i_walker_children_num > 1)
	{
		int i_came_from_child_index = -1;
		for (int ii = 0; ii < i_walker_children_num && i_came_from_child_index < 0; ii++)
		{
			if ((*pc_aux_walker)->vc_children[ii] == pc_came_from_child)
				i_came_from_child_index = ii;
		}// for (int ii = 0; ii < i_walker_children_num && i_came_from_child_index < 0; ii++)
		int i_rnd_way_to_go = rand() % (i_walker_children_num - 1);
		if (i_rnd_way_to_go == i_came_from_child_index)
			i_rnd_way_to_go = i_walker_children_num - 1;

		pc_aux_walker = &(*pc_aux_walker)->vc_children[i_rnd_way_to_go];
		bool bf_go_down = true;
		i_walker_children_num = (*pc_aux_walker)->vc_children.size();
		d_rnd = c_UNIF_REAL_PROB_DIST(c_GENERATOR);
		if (d_rnd < 1.0 - d_GO_UP_DOWN_PROB || i_walker_children_num == 0)
		{
			bf_go_down = false;
		}// if (d_rnd < 1.0 - d_GO_UP_DOWN_PROB || pc_aux_walker->pc_parent == nullptr)	
		while (bf_go_down)
		{
			i_walker_children_num = (*pc_aux_walker)->vc_children.size();
			i_rnd_way_to_go = rand() % i_walker_children_num;
			pc_aux_walker = &(*pc_aux_walker)->vc_children[i_rnd_way_to_go];
			i_walker_children_num = (*pc_aux_walker)->vc_children.size();
			d_rnd = c_UNIF_REAL_PROB_DIST(c_GENERATOR);
			if (d_rnd < 1.0 - d_GO_UP_DOWN_PROB || i_walker_children_num == 0)
			{
				bf_go_down = false;
			}// if (d_rnd < 1.0 - d_GO_UP_DOWN_PROB || pc_aux_walker->pc_parent == nullptr)	
		}// while (bf_go_down)
		//std::cout << "\npc_aux_walker: " << (*pc_aux_walker)->s_value << std::endl;  //debug
		CNode * pc_temp = pcSelected;
		pcSelected = *pc_aux_walker;
		*pc_aux_walker = pc_temp;
	}// if (i_walker_children_num > 1)
	else  // it means that we have reached root and root has just one child, just replace root with random generated tree (node)
	{
		if (*pc_aux_walker == pc_root)
			pc_aux_walker = &pc_root;
		//std::cout << "\npc_aux_walker: " << (*pc_aux_walker)->s_value << std::endl;  //debug
		//std::cout << "is it root?: " << (*pc_aux_walker == pc_root) << std::endl;  // debug
		v_replace_with_new_node(*pc_aux_walker);
	}// else
}// void CTree::v_swap_with_self_node(CNode * pcSelected)

CNode** CTree::pc_breadth_select(int iNodeNumber)
{
	std::queue<CNode **> q_nodes_queue;
	q_nodes_queue.push(&pc_root);
	int i_counter = 0;
	CNode** pc_rnd = nullptr;
	int i_node_children_num;
	CNode * pc_current_node;
	while (!q_nodes_queue.empty() && i_counter <= iNodeNumber)
	{
		if (i_counter == iNodeNumber)
		{
			pc_rnd = q_nodes_queue.front();
		}// if (i_counter == iNodeNumber)
		else
		{
			pc_current_node = *q_nodes_queue.front();
			i_node_children_num = pc_current_node->vc_children.size();
			for (int ii = 0; ii < i_node_children_num; ii++)
			{
				q_nodes_queue.push(&pc_current_node->vc_children[ii]);
			}// for (int ii = 0; ii < i_children_num; ii++)
		}// else
		q_nodes_queue.pop();
		i_counter++;
	}// while (!q_nodes_queue.empty() && i_counter <= iNodeNumber)
	return pc_rnd;
}// CNode * CTree::pc_breadth_select(int iNodeNumber)

void CTree::v_make_vars_and_leaves()
{
	vc_leaves.clear();
	msd_vars.clear();
	i_num_of_nodes = 0;
	v_var_leaves_making_walk(pc_root);
}// void CTree::v_make_vars_vec_and_leaves_vec()

void CTree::v_var_leaves_making_walk(CNode * pcNode)
{
	i_num_of_nodes++;
	int i_children_num = pcNode->vc_children.size();
	if (i_children_num == 0)
	{
		if (pcNode->b_is_var())
		{
			msd_vars.insert(std::pair<std::string, double>(pcNode->s_value, 0));
		}// if (pcNode->b_is_var())
		vc_leaves.push_back(pcNode);
	}// if (i_children_num == 0)
	for (int ii = 0; ii < i_children_num; ii++)
	{
		v_var_leaves_making_walk(pcNode->vc_children[ii]);
	}// for (int ii = 0; ii < i_children_num; ii++)
}// void CTree::v_var_leaves_vec_making_walk(CNode * pcNode)

void CTree::v_initialize_rnd_seed()
{
	if (!b_srand_was_called)
	{
		b_srand_was_called = true;
		srand(time(NULL));
	}// if (!b_srand_was_called)
}// void CTree::v_initialize_rnd_seed()

bool CTree::b_token_is_operator(std::string sToken)
{
	return (sToken == ADDITION || sToken == SUBTRACTION || sToken == MULTIPLICATION || sToken == DIVISION || sToken == SIN || sToken == COS);
}// bool CTree::b_token_is_operator(std::string sToken)

bool CTree::b_token_is_var(std::string sToken)
{
	int i_num_of_letters = 0;
	int i_num_of_numbers = 0;
	int i_token_len = sToken.length();
	bool b_is_var = true;
	if (i_token_len == 0)
		b_is_var = false;
	for (int ii = 0; ii < i_token_len && b_is_var; ii++)
	{
		if (sToken[ii] >= 48 && sToken[ii] <= 57)
		{
			i_num_of_numbers++;
		}// if (sToken[ii] >= 48 && sToken[ii] <= 57)
		else if ((sToken[ii] >= 65 && sToken[ii] <= 90) || (sToken[ii] >= 97 && sToken[ii] <= 122))
		{
			i_num_of_letters++;
		}// else if ((sToken[ii] >= 65 && sToken[ii] <= 90) || (sToken[ii] >= 97 && sToken[ii] <= 122))
		else
		{
			b_is_var = false;
		}// else
	}// for (int ii = 0; ii < i_token_len && b_is_var; ii++)
	if (b_is_var && i_num_of_letters == 0)
		b_is_var = false;

	return b_is_var;
}// bool CTree::b_token_is_var(std::string sToken)

std::string CTree::s_gen_rnd_operator()
{
	int i_operator_num = mis_CODE_TO_OPERATOR.size();
	int i_rnd = rand() % i_operator_num;
	return mis_CODE_TO_OPERATOR.at(i_rnd);
}// CNode * CTree::pc_gen_rnd_operator()

std::string CTree::s_gen_rnd_var_or_val()
{
	double d_rnd = c_UNIF_REAL_PROB_DIST(c_GENERATOR);
	std::string s_var_or_val;
	if (d_rnd >= d_VAR_PROBABILITY)
	{
		d_rnd = c_UNIF_REAL_PROB_DIST(c_GENERATOR);
		if (d_rnd >= d_X_VAR_PROBABILITY)  // 50% for X variable and 50 % for y variable
		{
			s_var_or_val = c_X_VAR;
		}// if (d_rnd >= 0.5)
		else
		{
			s_var_or_val = c_Y_VAR;
		}// else
	}// if (d_rnd >= 0.5)
	else  // value case
	{
		d_rnd = c_UNIF_REAL_PROB_DIST(c_GENERATOR);
		if (d_rnd >= 1.0 - d_DOUBLE_VAL_PROB)
		{
			double d_rnd_val = c_UNIF_DOUBLE_VAL_DIST(c_GENERATOR);
			s_var_or_val = std::to_string(d_rnd_val);
		}// if (d_rnd >= 1.0 - d_DOUBLE_VAL_PROB)
		else
		{
			int i_rnd = rand() % (i_RND_UPPER_BOUND_INC + 1 - i_RND_LOWER_BOUND_INC) + i_RND_LOWER_BOUND_INC;
			s_var_or_val = std::to_string(i_rnd);
		}// else
	}// else
	return s_var_or_val;
}// std::string CTree::s_gen_rnd_var_or_val()

std::string CTree::sNextToken(std::string const & sFormula, int *piCurrentPos)
{
	int i_substr_length = sFormula.length();
	int i_arg_end = 0;
	int i_arg_start = 0;
	bool blnf_arg_start = false;
	bool blnf_arg_end = false;
	for (int ii = *piCurrentPos; ii < i_substr_length && !blnf_arg_end; ii++)
	{
		if (!blnf_arg_start && sFormula[ii] != cSPACE)
		{
			i_arg_start = ii;
			blnf_arg_start = true;
		}// if (!blnf_arg_start && sFormula[ii] != ' ')
		if (blnf_arg_start && ((sFormula[ii + 1] == NULL) || (sFormula[ii + 1] == cSPACE)))
		{
			i_arg_end = ii;
			blnf_arg_end = true;
		}// if (blnf_arg_start && ((sFormula[ii + 1] == '\0') || (sFormula[ii + 1] == ' ')))
	}// for (int ii = 0; ii < i_substr_length && !blnf_arg_end; ii++)

	*piCurrentPos = i_arg_end + 1;
	return sFormula.substr(i_arg_start, i_arg_end - i_arg_start + 1);
}// std::string CTree::sNextToken(std::string const &sFormula, int *piCurrentPos)

int CTree::iNumOfTokensInString(std::string sFormula)
{
	int i_num_of_tokens = 0;
	int i_curr_pos = 0;
	int *pi_curr_pos = &i_curr_pos;
	int i_formula_len = sFormula.length();
	while (i_curr_pos < i_formula_len)
	{
		sNextToken(sFormula, pi_curr_pos);
		i_num_of_tokens++;
	}// while (i_curr_pos < i_formula_len)
	return i_num_of_tokens;
}// int CTree::iNumOfTokensInString(std::string sFormula)

bool CTree::b_srand_was_called = false;


// CNode part


CNode::CNode()
{
	pc_parent = nullptr;
	i_depth = 0;
}// CNode::CNode()

CNode::CNode(int iMaxHeight, CNode* pcParent)
{
	pc_parent = pcParent;
	i_depth = pcParent->i_depth + 1;
}// CNode::CNode(int iParentDepth, CNode* pcParent)

CNode::CNode(CNode * pcParent)
{
	this->pc_parent = pcParent;
}// CNode::CNode(CNode * pcParent)

CNode::~CNode()
{
	int i_children_number = vc_children.size();
	for (int ii = 0; ii < i_children_number; ii++)
	{
		delete vc_children[ii];
	}// for (int ii = 0; ii < i_children_number; ii++)
	vc_children.clear();
}// CNode::~CNode()

bool CNode::b_is_var()
{
	bool b_is_var = false;
	if (vc_children.empty() && !CTree::bIsDouble(s_value) && !CTree::bIsInteger(s_value))
	{
		b_is_var = true;
	}// if (vc_children.empty() && !bIsInteger(s_value))
	return b_is_var;
}// bool CNode::b_is_var()


// functions

bool CTree::bIsInteger(std::string sToCheck)
{
	int i_str_length = sToCheck.length();
	bool blnf_is_number = true;
	if (sToCheck[0] == 48 && i_str_length > 1)
		blnf_is_number = false;
	if (sToCheck[0] < 48 || sToCheck[0] > 57)
		blnf_is_number = false;
	if (sToCheck[0] == 45)
		blnf_is_number = true;
	for (int ii = 1; ii < i_str_length && blnf_is_number; ii++)
	{
		if (sToCheck[ii] < 48 || sToCheck[ii] > 57)
			blnf_is_number = false;
	}//for (int ii = 0; ii < i_str_length && blnf_is_number; ii++)
	if ((sToCheck[0] == 45 && i_str_length == 1) || sToCheck[0] == NULL || sToCheck == sEMPTY)
		blnf_is_number = false;
	return blnf_is_number;
}// bool bIsInteger(std::string const sToCheck)

bool CTree::bIsDouble(std::string sToCheck)
{
	int i_str_length = sToCheck.length();
	bool blnf_is_number = true;
	if (i_str_length > 0)
	{
		if (sToCheck[0] == 48 && i_str_length > 1 && sToCheck[1]!=cDOT && sToCheck[1]!=cCOMMA)
			blnf_is_number = false;
		if (sToCheck[0] < 48 || sToCheck[0] > 57)
			blnf_is_number = false;
		if (sToCheck[0] == 45)
			blnf_is_number = true;
		if (sToCheck[0] == 45 && i_str_length > 1 && (sToCheck[1] < 48 || sToCheck[1]>57))
			blnf_is_number = false;
		if (sToCheck[0] == 45 && i_str_length > 2 && sToCheck[1] == 48 && sToCheck[2] == 48)
			blnf_is_number = false;
		bool blnf_dot_appeard = false;
		for (int ii = 1; ii < i_str_length && blnf_is_number; ii++)
		{
			if ((sToCheck[ii] == cDOT || sToCheck[ii] == cCOMMA) && !blnf_dot_appeard)
				blnf_dot_appeard = true;
			else if (sToCheck[ii] < 48 || sToCheck[ii] > 57)
				blnf_is_number = false;
		}//for (int ii = 0; ii < i_str_length && blnf_is_number; ii++)
		if (blnf_is_number && (sToCheck[0] == 45 && i_str_length == 1 || sToCheck[i_str_length - 1] == cDOT || sToCheck[i_str_length - 1] == cCOMMA))
			blnf_is_number = false;
		if (!blnf_is_number && i_str_length == 2 && sToCheck[0] == 43 && sToCheck[1] == 48)
			blnf_is_number = true;
		for (int ii = 0; ii < i_NUM_OF_UNUSUAL_DOUBLES && !blnf_is_number; ii++)
		{
			if(sToCheck == pc_POSSIBLE_UNUSUAL_DOUBLES[ii])
				blnf_is_number = true;
		}// for (int ii = 0; ii < i_NUM_OF_UNUSUAL_DOUBLES && !blnf_is_number; ii++)	
	}// if (i_str_length > 0)
	else
	{
		blnf_is_number = false;
	}// else
	return blnf_is_number;
}// bool bIsDouble(std::string sToCheck)

