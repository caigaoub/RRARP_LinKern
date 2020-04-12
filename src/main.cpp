#include <stdlib.h>
#include "DataHandler.h"
#include "PartitionScheme.h"
#include <ctime>
#include <chrono>
#include <boost/lexical_cast.hpp>
#include <sys/stat.h>
#include "PulseAlgo.h"
#include "ProgTime.h"
using namespace std;

int main(int argc, const char* argv[]) {
	argc = argc; // get rid of warning: unused argc
	// const int nb_obps = atoi(argv[1]);
	// const int type_trajc = atoi(argv[2]);
	string configfile = argv[1];
	const int nb_obps = atoi(argv[2]);
	const int source = atoi(argv[3]);
	const int sink = atoi(argv[4]);

	const int log_on = atoi(argv[5]);

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
	// dataset_.print();
	Pulse  test_;
	// test_.initialize(dataset_, 4, 16, nb_obps);
	test_.initialize(dataset_, 4, 16, nb_obps, cur_dir+"InnerGraphs/n_6_h_1.graph_"+to_string(nb_obps), source, sink);

	// Vertex entry = test_.get_turnpoint(0);
	// Vertex exit = test_.get_turnpoint(2);
	// test_.update_graph(entry, exit);
	test_.calc_leastrisk_sink();
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
	vector<int> path;
	ProgTime test_time;
	test_time.start_prog();
	test_.recursive_search(0, path, 0, 0, log_on);
	test_time.end_prog();
	cout << " =====>>>> Time: " << test_time._elapsed_secs << endl;
	test_.print_opt_sol();
//	system("pause");
	return 0;
}



