#include <stdlib.h>
#include "DataHandler.h"
#include <ctime>
#include <chrono>
#include <boost/lexical_cast.hpp>
#include <sys/stat.h>
#include <vector>
#include <thread>
#include <fstream>
#include <sstream>
#include "LKAlgo.h"
#include "TSPModels.h"

using namespace std;


int main(int argc, const char* argv[]) {
	argc = argc; // get rid of warning: unused argc
	string filename = argv[1];
	const int nb_targets = atoi(argv[2]);


	DataHandler dataset_;
	// string filename = "/home/cai/Dropbox/Box_Research/Github/RRARP_LinKern/dat/Inst_LK/"
	dataset_.read_instance(filename);
	dataset_.build_riskgraph();


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



	TSPModel_STE tsp_model_;
	tsp_model_.init_edge_weights(nb_targets+2, dataset_._c2c_adjmatr);
	tsp_model_.create_formula();
	tsp_model_.solve();
	auto tspsol = test_.solve_shortestpath_path(tsp_model_._opt_seq);


	TSPModel_STE lb_model_;
	lb_model_.init_edge_weights(nb_targets+2, dataset_._b2b_adjmatr);
	lb_model_.create_formula();
	lb_model_.solve();

	cout << "LB: " << dataset_._total_inrisk_LB + lb_model_._curbestobj << endl;
	cout << "tsp sol: " << tspsol.second << endl;
	cout << "LK sol: " << test_._Obj_Star << endl;

	return 0;
}


int main2(int argc, const char* argv[]) {
	argc = argc; // get rid of warning: unused argc
	string filename = argv[1];
	const int nb_targets = atoi(argv[2]);


	DataHandler dataset_;
	// string filename = "/home/cai/Dropbox/Box_Research/Github/RRARP_LinKern/dat/Inst_LK/"
	dataset_.read_instance(filename);
	dataset_.build_riskgraph();




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


