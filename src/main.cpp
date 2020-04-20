#include <stdlib.h>
#include "DataHandler.h"
#include "PartitionScheme.h"
#include <ctime>
#include <chrono>
#include <boost/lexical_cast.hpp>
#include <sys/stat.h>
#include <vector>
#include "PulseAlgo.h"
#include "ProgTime.h"
#include<thread>
using namespace std;

void func_parr(DataHandler& dataset_, int nb_obps, string cur_dir, int source, int sink, double demandPCT){
	// cout << "-------------------------------------" << endl;
	// cout << "i = " << source << ", j = " << sink << endl;
	Pulse  test_;

	test_.initialize_triGraph(dataset_, 4, nb_obps, cur_dir+"InnerGraphs/n_6_h_1.graph_"+to_string(nb_obps));
	test_.set_parameters(source, sink, demandPCT);
	test_.calc_leastrisk_sink();
	vector<int> path;
	// ProgTime test_time;
	// test_time.start_prog();
	test_.recursive_search(source, path, 0, 0, 0);	
	// test_time.end_prog();
	// cout << " =====>>>> Time: " << test_time._elapsed_secs << endl;
	// test_.print_opt_sol();
}


int main(int argc, const char* argv[]) {
	argc = argc; // get rid of warning: unused argc
	// const int nb_obps = atoi(argv[1]);
	// const int type_trajc = atoi(argv[2]);
	string configfile = argv[1];
	const int nb_obps = atoi(argv[2]);
	// const int source = atoi(argv[3]);
	// const int sink = atoi(argv[4]);
	const double demandPCT = atof(argv[3]);
	// const int log_on = atoi(argv[4]);


	fstream file(configfile);
	if (!file) {
		cerr << "ERROR: could not open config '" << configfile << "' for reading'" << endl;
		throw(-1);
	}
	string instance_name_only;
	file >> instance_name_only;
	file.close();
	// string cur_dir = "/projects/academic/josewalt/caigao/RRARP_LinKern/dat/";
    string cur_dir = "/home/cai/Dropbox/Box_Research/Github/RRARP_LinKern/dat/";
	struct stat buffer;
  	if(stat (cur_dir.c_str(), &buffer) != 0){
  		cerr << " Path of instances does not exist!! (in main.cpp:line 29) " << endl;
  	}
	string instance_wPath = cur_dir + instance_name_only;
	DataHandler dataset_;
	dataset_.parse(instance_wPath);

    /*===================== start parallel matching  ================================**/
    unsigned n_available = std::thread::hardware_concurrency();
    int nr_threads = (int)(n_available-1);
//    unsigned nr_threads = 1;
    vector<pair<int,int>> idxset;
    for(int i =0; i<16; i++){
    	for(int j=0; j<16; j++){
    		if(i != j)
    			idxset.push_back(make_pair(i,j));
    	}
    }
    for(int i=0; i < idxset.size(); i++)
    	cout << idxset[i].first << ", " << idxset[i].second << endl;
    int itr =0;
    int size = (int)idxset.size();
    while(itr < size){
    	vector<thread> threads;
    	int cores = min(nr_threads, size - itr);
		for (int k = 0; k < cores; ++k) {
    	  	cout << idxset[itr].first << ", " << idxset[itr].second <<endl;
		    threads.push_back(thread(func_parr, ref(dataset_), nb_obps, cur_dir, idxset[itr].first, idxset[itr].second, demandPCT));
		    itr++;
	    }
   	    for(auto &t : threads){ // Join the threads with the main thread
	        t.join();
	    }
    }



	// for(unsigned i=0; i<test_._lbrisk_sink.size(); i++)
	// 	cout << test_._lbrisk_sink[i] << ' ';
	// vector<vector<double>> adjmatx({{0, 2, 3, 0, 0, 0, 0}, 
	// 					 {2, 0, 9, 1, 2, 0, 0}, 
	// 					 {3, 9, 0, 1, 0, 1, 0}, 
	// 					 {0, 1, 1, 0, 1, 1, 0}, 
	// 					 {0, 2, 0, 1, 0, 5, 2},
	// 					 {0, 0, 1, 1, 5, 0, 1},
	// 					 {0, 0, 0, 0, 2, 1, 0}});
	// int size = 7;
	// vector<double> rewards({0, 2, 3, 3, 4, 2, 0});
	// test_.initialize_generalgraph(adjmatx, rewards, size, 0.5);
	// test_.calc_leastrisk_sink();


	// test_.set_parameters(source, sink, demandPCT);
	// test_.calc_leastrisk_sink();

	// vector<int> path;
	// ProgTime test_time;
	// test_time.start_prog();
	// test_.recursive_search(source, path, 0, 0, log_on);
	// test_time.end_prog();
	// cout << " =====>>>> Time: " << test_time._elapsed_secs << endl;
	// test_.print_opt_sol();
//	system("pause");
	return 0;
}



