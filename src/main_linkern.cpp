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

using namespace std;


int main(int argc, const char* argv[]) {
	// argc = argc; // get rid of warning: unused argc
	// // const int nb_obps = atoi(argv[1]);
	// // const int type_trajc = atoi(argv[2]);
	// string configfile = argv[1];
	// string configfile2 = argv[2];

	// // const int nb_obps = atoi(argv[2]);
	// // const int source = atoi(argv[3]);
	// // const int sink = atoi(argv[4]);
	// // const double demandPCT = atof(argv[3]);
	// const int is_server = atoi(argv[3]);

	// fstream file_params(configfile2);
	// if (!file_params) {
	// 	cerr << "ERROR: could not open config '" << configfile2 << "' for reading'" << endl;
	// 	throw(-1);
	// }
	// int nb_obps; 
	// double demandPCT;
	// file_params >> nb_obps >> demandPCT;
	// file_params.close();

	// fstream file(configfile);
	// if (!file) {
	// 	cerr << "ERROR: could not open config '" << configfile << "' for reading'" << endl;
	// 	throw(-1);
	// }
	// string instance_name_only;
	// file >> instance_name_only;
	// file.close();
	// string cur_dir = " ";
	// if(is_server){
	// 	cur_dir = "/projects/academic/josewalt/caigao/RRARP_LinKern/dat/";
	// }else{
	//     cur_dir = "/home/cai/Dropbox/Box_Research/Github/RRARP_LinKern/dat/";
	// }
	// struct stat buffer;
 //  	if(stat (cur_dir.c_str(), &buffer) != 0){
 //  		cerr << " Path of instances does not exist!! (in main.cpp:line 29) " << endl;
 //  	}
	// string instance_wPath = cur_dir + instance_name_only;
	// DataHandler dataset_;
	// dataset_.parse(instance_wPath);


 //    vector<pair<int,int>> idxset;
 //    for(int i =0; i<16; i++){
 //    	for(int j=i+1; j<16; j++){
 //    		idxset.push_back(make_pair(i,j));
 //    	}
 //    }
 //    double total_time = 0.0;
 //    for(unsigned int i=0; i<idxset.size(); i++){
 //    	cout << idxset[i].first << ", " << idxset[i].second << endl;
	// 	double temp_time = func_parr(dataset_, nb_obps, cur_dir, idxset[i].first, idxset[i].second, demandPCT);
	// 	total_time += temp_time;
 //    }

 //    total_time /= (double)idxset.size();

 //    cout << "time: " << total_time <<  endl;

	// ofstream file2;
	// string out_dir = " ";
	// if(is_server){
	// 	out_dir = "/projects/academic/josewalt/caigao/RRARP_LinKern/ret/pulse_out/";
	// }else{
	// 	out_dir	= "/home/cai/Dropbox/Box_Research/Github/RRARP_LinKern/ret/pulse_out/";
	// }
 //    file2.open(out_dir + "triGraph_" + to_string(nb_obps) + "_" +to_string(demandPCT));
	// file2 << "time: " << total_time << '\n';
	// file2.close();


//     /*===================== start parallel matching  ================================**/
//     unsigned n_available = std::thread::hardware_concurrency();
//     int nr_threads = (int)(n_available-1);
// //    unsigned nr_threads = 1;
//     vector<pair<int,int>> idxset;
//     for(int i =0; i<16; i++){
//     	for(int j=0; j<16; j++){
//     		if(i != j)
//     			idxset.push_back(make_pair(i,j));
//     	}
//     }
//     for(int i=0; i < idxset.size(); i++)
//     	cout << idxset[i].first << ", " << idxset[i].second << endl;
//     int itr =0;
//     int size = (int)idxset.size();
//     while(itr < size){
//     	vector<thread> threads;
//     	int cores = min(nr_threads, size - itr);
// 		for (int k = 0; k < cores; ++k) {
//     	  	cout << idxset[itr].first << ", " << idxset[itr].second <<endl;
// 		    threads.push_back(thread(func_parr, ref(dataset_), nb_obps, cur_dir, idxset[itr].first, idxset[itr].second, demandPCT));
// 		    itr++;
// 	    }
//    	    for(auto &t : threads){ // Join the threads with the main thread
// 	        t.join();
// 	    }
//     }




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



