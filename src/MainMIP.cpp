#include <stdlib.h>
#include "gurobi_c++.h"
#include "DataHandler.h"
#include "STEFormulation.h"
#include "SubtourCuts.h"
#include <ctime>
#include <chrono>
#include <boost/lexical_cast.hpp>
#include <sys/stat.h>
#include <vector>
#include <thread>
#include <fstream>
#include <sstream>
#include "LKAlgo.h"
using namespace std;

int main(int argc, const char* argv[]) {
	argc = argc; // get rid of warning: unused argc
	string filename = argv[1];
	const int nb_targets = atoi(argv[2]);


	DataHandler dataset_;
	// string filename = "/home/cai/Dropbox/Box_Research/Github/RRARP_LinKern/dat/Inst_LK/"
	dataset_.read_instance(filename);
	dataset_.build_riskgraph();
	/*Gurobi model for master problem */
	GRBEnv * evn_MP_ = new GRBEnv();
	GRBModel model_MP_ = GRBModel(*evn_MP_);
	STEFormulation formul_master;
	// model_MP_.getEnv().set(GRB_IntParam_OutputFlag, 0);
	formul_master.build_formul(&model_MP_, dataset_);
	formul_master.solve_formul_wCB();


	exit(0);		

	LK test_;
	test_.initialize(dataset_);

	// vector<int> seq = {0,1,2,3,4,5,6,7};
	// vector<int> seq = {0,5,1,3,2,4,6,7};
	// vector<int> seq = {0,4,1,5,3,2,6,7};
	// vector<int> seq = {0,6,2,1,3,4,5,7};
	// vector<int> seq1 = {0,1,2,3,4,5,6,7,8,9,10,11};
	// vector<int> seq1 = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21};
	vector<int> seq1;
	seq1.push_back(0);
	for(int i=0; i<nb_targets; i++)
		seq1.push_back(i+1);
	seq1.push_back(nb_targets+1);

	// vector<int> seq1 = {0, 9, 10, 5, 6, 7, 8, 3, 2, 1, 4, 11};
	
	// vector<int> seq1 = {0, 10, 9, 4, 8, 2, 3, 7, 1, 6, 5, 11};
	auto initialtour1 = test_.solve_shortestpath_path(seq1);
	test_.solve(initialtour1.first, false);

	


	return 0;
}