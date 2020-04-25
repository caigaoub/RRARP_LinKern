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
using namespace std;


int main(int argc, const char* argv[]) {
	argc = argc; // get rid of warning: unused argc
	string filename = argv[1];

	DataHandler dataset_;
	// string filename = "/home/cai/Dropbox/Box_Research/Github/RRARP_LinKern/dat/Inst_LK/"
	dataset_.read_instance(filename);
	dataset_.build_riskgraph();


	LK test_;
	test_.initialize(dataset_);

	// vector<int> seq = {0,1,2,3,4,5,6,7};
	// vector<int> seq = {0,5,1,3,2,4,6,7};
	vector<int> seq = {0,4,1,5,3,2,6,7};
	TOUR initialtour = test_.solve_shortestpath_path(seq);
	// test_.random_tsp_seq(16);
	test_.solve(initialtour, true);

	return 0;
}



