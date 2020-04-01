#include <stdlib.h>
#include "DataHandler.h"
#include "PartitionScheme.h"
#include <ctime>
#include <chrono>
#include <boost/lexical_cast.hpp>
#include <sys/stat.h>


using namespace std;

int main(int argc, const char* argv[]) {
	argc = argc; // get rid of warning: unused argc
	const int nb_dstzn = atoi(argv[1]);
	const int type_trajc = atoi(argv[2]);
	string configfile = argv[3];

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
	dataset_.print();
	PartitionScheme network_;
	network_.build(dataset_, nb_dstzn, type_trajc);

//	system("pause");
	return 0;
}

