#pragma once
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <ctime>
#include <random>
#include <queue>

#define cFIX_NUM '1'
#define cSPACE ' '
#define sEMPTY ""
#define cDOT '.'
#define cCOMMA ','
#define sUNALLOWED_EXP "unallowed expression "
#define ADDITION "+"
#define SUBTRACTION "-"
#define MULTIPLICATION "*"
#define DIVISION "/"
#define SIN "sin"
#define COS "cos"

static const int i_NUM_OF_UNUSUAL_DOUBLES = 8;
static const std::string pc_POSSIBLE_UNUSUAL_DOUBLES [i_NUM_OF_UNUSUAL_DOUBLES] = {"nan", "-nan", "-nan(ind)", "nan(ind)", "inf", "-inf", "+0", "-0"};
static std::default_random_engine c_GENERATOR(time(NULL));
static const std::uniform_real_distribution<double> c_UNIF_REAL_PROB_DIST;
static const std::map<int, std::string> misCodeToOperInit();
static const std::map<int, std::string> mis_CODE_TO_OPERATOR = misCodeToOperInit();
static const std::map<std::string, int> msiOpToChildrenNumInit();
static const std::map<std::string, int> msi_OPERATOR_CHILDREN_NUM = msiOpToChildrenNumInit();
static const double d_OPERATOR_GENERATE_PROBABILITY = 0.4;
static const double d_ROOT_OPERATOR_PROBABILITY = 0.9;
static const double d_VAR_PROBABILITY = 0.6;
static const double d_X_VAR_PROBABILITY = 0.5;
static const double d_DOUBLE_VAL_PROB = 0.75;
static const char c_X_VAR = 'x';
static const char c_Y_VAR = 'y';
static const int i_RND_LOWER_BOUND_INC = -10;
static const int i_RND_UPPER_BOUND_INC = 10;
static const std::uniform_real_distribution<double> c_UNIF_DOUBLE_VAL_DIST(i_RND_LOWER_BOUND_INC, i_RND_UPPER_BOUND_INC);
static const double d_SELF_NODE_SWAP_PROB = 0.3;
static const double d_GO_UP_DOWN_PROB = 0.65;
static const int i_MAX_HEIGHT = 3;

class CNode;
enum class CTreeError;

class CTree
{
	friend int main();  // testing purpose
public:
	CTree(CTree const & cOtherTree);
	CTree(int iMaxHeight);
	~CTree();
	static CTree * pcConstructFromString(std::string sFormula);
	std::string sGetFormula();
	std::string sGetVars();
	std::string sGetLeaves();  // testing purpose
	double dEvaluateFormula(std::vector<double> const &c_values);
	CTree operator+(CTree const &cOtherTree);
	CTree & operator=(CTree const &cOtherTree);
	CTreeError errGetErrorCode();
	void vMutation();
	static void vCrossover(CTree & cTree1, CTree & cTree2);
	double dGetDistance();
	void vSetDistance(double dDistance);
	bool bGetChangeStatus();
	void vSetIfChanged(bool b_changed);
	int iGetNumOfNodes();
	static bool bGetSrandCallingStatus();
	static void vSetSrandCallingStatus(bool bStatus);
	static bool bIsInteger(std::string sToCheck);
	static bool bIsDouble(std::string sToCheck);

private:
	CTree();
	CTree(std::string sInputFormula);
	void v_construct_node(CNode * pcNode, std::string const & sFormula, int *piCurrentPos);
	void v_copy_node(CNode * pcNode, CNode & cOtherNode);
	static std::string s_fix_input_formula(std::string sFormula);
	std::string s_aux_preorder_walk(CNode *pcNode);
	double d_aux_evaluate(CNode *pcNode);
	double d_get_value_of_var(std::string const & sVarName);
	CNode* pc_gen_rnd_root(int iMaxHeight);
	
	void v_gen_rnd_tree(int iMaxHeight, int iRootChildrenNum);
	void v_gen_rnd_node(int iMaxHeight, CNode * pcNode);
	CNode** pc_get_rnd_node(bool bfRootIsAcceptable);
	void v_replace_with_new_node(CNode*& pcSelected);
	void v_swap_with_self_node(CNode*& pcSelected);
	CNode** pc_breadth_select(int iNodeNumber);
	void v_make_vars_and_leaves();
	void v_var_leaves_making_walk(CNode *pcNode);
	void v_initialize_rnd_seed();

	static bool b_token_is_operator(std::string sToken);
	static bool b_token_is_var(std::string sToken);
	static std::string s_gen_rnd_operator();
	static std::string s_gen_rnd_var_or_val();
	static std::string sNextToken(std::string const & sFormula, int *piCurrentPos);
	static int iNumOfTokensInString(std::string sFormula);

	CNode *pc_root;
	std::map<std::string, double> msd_vars;
	std::vector<CNode *> vc_leaves;
	CTreeError err_code;
	int i_num_of_nodes;
	double d_distance_measure;
	bool b_changed_since_eval;

	static bool b_srand_was_called;
};// class CTree

class CNode
{
	// all private, this class is only for CTree usage
	friend class CTree;

	CNode();
	CNode(int iMaxHeight, CNode *pcParent);
	CNode(CNode * pcParent);
	~CNode();

	bool b_is_var();

	std::string s_value;
	CNode *pc_parent;
	std::vector<CNode*> vc_children;
	int i_depth;
};// class CNode

enum class CTreeError
{
	ERR_OK = 0,
	ERR_TO_MANY_VALUES = 1,
	ERR_NOT_ENOUGH_VALUES = 2,
};// enum class CTreeError

const std::map<int, std::string> misCodeToOperInit()
{
	std::map<int, std::string> mis_code_top_op;
	mis_code_top_op[0] = ADDITION;
	mis_code_top_op[1] = SUBTRACTION;
	mis_code_top_op[2] = MULTIPLICATION;
	mis_code_top_op[3] = DIVISION;
	mis_code_top_op[4] = SIN;
	mis_code_top_op[5] = COS;
	return mis_code_top_op;
}// const std::map<int, std::string> misCodeToOperInit()

const std::map<std::string, int> msiOpToChildrenNumInit()
{
	std::map<std::string, int> msi_op_to_children_num;
	msi_op_to_children_num[ADDITION] = 2;
	msi_op_to_children_num[SUBTRACTION] = 2;
	msi_op_to_children_num[MULTIPLICATION] = 2;
	msi_op_to_children_num[DIVISION] = 2;
	msi_op_to_children_num[SIN] = 1;
	msi_op_to_children_num[COS] = 1;
	return msi_op_to_children_num;
}// const std::map<std::string, int> msiOpToChildrenNumInit()