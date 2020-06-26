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
#include <numeric>
using namespace std;
#include <iomanip>
void print(vector<double>& seq){
	// cout << "tsp: ";
	for(auto itr=seq.begin(); itr != seq.end(); itr++){
		cout << *itr << " & ";
	}
	cout << '\n';
}
int main(int argc, const char* argv[]) {
	argc = argc; // get rid of warning: unused argc
	string filename = argv[1];
	const int nb_targets = atoi(argv[2]);
	DataHandler dataset_;
	// string filename = "/home/cai/Dropbox/Box_Research/Github/RRARP_LinKern/dat/Inst_LK/"
	dataset_.read_instance(filename);

	dataset_.build_riskgraph();


	LK testLK_;
	testLK_.initialize(dataset_);
	LK testOpt2_;
	testOpt2_.initialize(dataset_);

	int maxround = 20;
	vector<double> objsLK;
	vector<double> objsOpt2;
	vector<double> timeLK;
	vector<double> timeOpt2;
	for(int i =0; i < maxround; i++){
		auto seq = testLK_.random_tsp_seq(nb_targets);
		testLK_.print(seq);
		auto initialtour1 = testLK_.solve_shortestpath_path(seq);
		ProgTime lk;
		lk.start_prog();
		testLK_.solve_linkern(initialtour1.first, false);
		lk.end_prog();
		timeLK.push_back(lk._elapsed_secs);
		cout << "LK time: " << lk._elapsed_secs << " ** " << testLK_._Obj_Star << endl;
		objsLK.push_back(testLK_._Obj_Star);


		// initialtour1 = testOpt2_.solve_shortestpath_path(seq);
		ProgTime opt2;
		opt2.start_prog();
		testOpt2_.solve_opt2(initialtour1.first, false);
		opt2.end_prog();
		timeOpt2.push_back(opt2._elapsed_secs);


		cout << "opt2 time: " << opt2._elapsed_secs << " ** " << testOpt2_._Obj_Star << endl;

		objsOpt2.push_back(testOpt2_._Obj_Star);
		exit(0);
	}

	// int maxround = 5;
	// set<double> objsOpt2;
	// for(int i =0; i < maxround; i++){
	// 	auto seq = testOpt2_.random_tsp_seq(nb_targets);
	// 	testOpt2_.print(seq);
	// 	auto initialtour1 = testOpt2_.solve_shortestpath_path(seq);
	// 	testOpt2_.solve_opt2(initialtour1.first, false);
	// 	objsOpt2.insert(testOpt2_._Obj_Star);
	// }
	print(objsLK);
	print(objsOpt2);

	double avgtime_lk = std::accumulate(timeLK.begin(), timeLK.end(), 0.0)/(double)maxround;
	double avgtime_opt2 = std::accumulate(timeOpt2.begin(), timeOpt2.end(), 0.0)/(double)maxround;
	


	// // cout << "LK sol: " << test_._Obj_Star << endl;


	TSPModel_STE tsp_model_;
	tsp_model_.init_edge_weights(nb_targets+2, dataset_._c2c_adjmatr);
	tsp_model_.create_formula();
	tsp_model_.solve();
	auto tspsol = testLK_.solve_shortestpath_path(tsp_model_._opt_seq);
	// string dir = ROOTDIR;
	// string filename2 = dir+"/dat/pyFunc/optimalpath_" + to_string(0) + ".txt";
	// test_.write_optimalpath(filename2, tspsol.first);

	TSPModel_STE lb_model_;
	lb_model_.init_edge_weights(nb_targets+2, dataset_._b2b_adjmatr);
	lb_model_.create_formula();
	lb_model_.solve();
	double LB = dataset_._total_inrisk_LB + lb_model_._curbestobj;
	double avgGapOPT2 = 0;
	for(auto itr=objsOpt2.begin(); itr != objsOpt2.end(); itr++ ){
		avgGapOPT2 += (*itr - LB)/LB * 100.0;
	}
	avgGapOPT2 /= (double)maxround;
	double avgGapLK = 0;
	for(auto itr=objsLK.begin(); itr != objsLK.end(); itr++ ){
		avgGapLK += (*itr - LB)/LB * 100.0;
	}
	avgGapLK /= (double)maxround;


	cout << std::fixed << std::setprecision(2) << LB << " & ";
	cout << tspsol.second << " & " << (tspsol.second- LB)/LB * 100.0 << " & ";
	cout << *min_element(objsOpt2.begin(), objsOpt2.end()) << " & " << *max_element(objsOpt2.begin(), objsOpt2.end())  << " & ";
	cout << avgGapOPT2 << " & " << avgtime_opt2  << " & ";
	cout << *min_element(objsLK.begin(), objsLK.end()) << " & " << *max_element(objsLK.begin(), objsLK.end())  << " & ";
	cout << avgGapLK << " & " << avgtime_lk << endl;

	// cout << "LK sol: " << test_._Obj_Star << endl;
	// testLK_.print(objsLK);
	// testLK_.print(objsOpt2);
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


	LK testLK_;
	testLK_.initialize(dataset_);
	LK testOpt2_;
	testOpt2_.initialize(dataset_);
	// vector<int> seq = {0,1,2,3,4,5,6,7};
	// vector<int> seq = {0,5,1,3,2,4,6,7};
	// vector<int> seq = {0,4,1,5,3,2,6,7};
	// vector<int> seq = {0,6,2,1,3,4,5,7};
	// vector<int> seq1 = {0,1,2,3,4,5,6,7,8,9,10,11};
	// vector<int> seq1 = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21};
	// vector<int> seq1;
	// seq1.push_back(0);
	// for(int i=0; i<nb_targets; i++)
	// 	seq1.push_back(i+1);
	// seq1.push_back(nb_targets+1);
	

	
	int maxround = 10;
	vector<double> objsLK;
	vector<double> objsOpt2;

	for(int i =0; i < maxround; i++){
		auto seq = testLK_.random_tsp_seq(nb_targets);
		testLK_.print(seq);
		auto initialtour1 = testLK_.solve_shortestpath_path(seq);
		ProgTime lk;
		lk.start_prog();
		testLK_.solve_linkern(initialtour1.first, false);
		lk.end_prog();
		cout << "LK time: " << lk._elapsed_secs << endl;
		objsLK.push_back(testLK_._Obj_Star);


		// initialtour1 = testOpt2_.solve_shortestpath_path(seq);
		ProgTime opt2;
		opt2.start_prog();
		testOpt2_.solve_opt2(initialtour1.first, false);
		opt2.end_prog();
		cout << "opt2 time: " << opt2._elapsed_secs << endl;

		objsOpt2.push_back(testOpt2_._Obj_Star);
	}

	// int maxround = 5;
	// set<double> objsOpt2;
	// for(int i =0; i < maxround; i++){
	// 	auto seq = testOpt2_.random_tsp_seq(nb_targets);
	// 	testOpt2_.print(seq);
	// 	auto initialtour1 = testOpt2_.solve_shortestpath_path(seq);
	// 	testOpt2_.solve_opt2(initialtour1.first, false);
	// 	objsOpt2.insert(testOpt2_._Obj_Star);
	// }
	print(objsLK);
	print(objsOpt2);

	double sum_lk = std::accumulate(objsLK.begin(), objsLK.end(), 0.0);
	double mean_lk = sum_lk / (double)objsLK.size();

	double sq_sum_lk = std::inner_product(objsLK.begin(), objsLK.end(), objsLK.begin(), 0.0);
	double stdev_lk = std::sqrt(sq_sum_lk / objsLK.size() - mean_lk * mean_lk);
	cout << mean_lk << " & " << stdev_lk << endl; 


	double sum_opt2 = std::accumulate(objsOpt2.begin(), objsOpt2.end(), 0.0);
	double mean_opt2 = sum_opt2 / (double)objsOpt2.size();

	double sq_sum_opt2 = std::inner_product(objsOpt2.begin(), objsOpt2.end(), objsOpt2.begin(), 0.0);
	double stdev_opt2 = std::sqrt(sq_sum_opt2 / objsOpt2.size() - mean_opt2 * mean_opt2);
	cout << mean_opt2 << " & " << stdev_opt2 << endl; 


	// // cout << "LK sol: " << test_._Obj_Star << endl;


	// TSPModel_STE tsp_model_;
	// tsp_model_.init_edge_weights(nb_targets+2, dataset_._c2c_adjmatr);
	// tsp_model_.create_formula();
	// tsp_model_.solve();
	// auto tspsol = testLK_.solve_shortestpath_path(tsp_model_._opt_seq);
	// string dir = ROOTDIR;
	// string filename2 = dir+"/dat/pyFunc/optimalpath_" + to_string(0) + ".txt";
	// // test_.write_optimalpath(filename2, tspsol.first);

	// TSPModel_STE lb_model_;
	// lb_model_.init_edge_weights(nb_targets+2, dataset_._b2b_adjmatr);
	// lb_model_.create_formula();
	// lb_model_.solve();

	// cout << "LB: " << dataset_._total_inrisk_LB + lb_model_._curbestobj << endl;
	// cout << "TSP sol: " << tspsol.second << endl;
	// // cout << "LK sol: " << test_._Obj_Star << endl;
	// testLK_.print(objsLK);
	// testLK_.print(objsOpt2);
	return 0;
}

