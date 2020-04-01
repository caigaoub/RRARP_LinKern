#include "DataHandler.h"


void DataHandler::parse(string filename) {
	fstream file(filename);
	if (!file) {
		cerr << "ERROR: could not open file '" << filename << "' for reading'" << endl;
		throw(-1);
	}
	// Geographic Information of Targets
	auto pos = filename.find_last_of("/");
    _name = filename.substr(pos+1, filename.size());
    size_t found = _name.find("c");
	if (found == string::npos) {
		file >> _nb_targets;
		cout << "This instance is clustered!" << endl;
	}else {
		file >> _nb_targets >> _nb_clusters;
		cout << "This is a clustered instance:" << _nb_clusters << endl;
	}
	file >> _depot1_loc._x >> _depot1_loc._y;
	file >> _depot2_loc._x >> _depot2_loc._y;

	_target_locs = new Vertex[_nb_targets];
	_radii = new double[_nb_targets];
	for (int i = 0; i < _nb_targets; i++) {
		file >> _target_locs[i]._x >> _target_locs[i]._y >> _radii[i];
	}

	_bdg_rewards_ratio = new double[_nb_targets];
	// _risk_thold_ratio = new double[_nb_targets];
	for (int i = 0; i < _nb_targets; i++) {
		file >> _bdg_rewards_ratio[i];
	}
	// for (int i = 0; i < _nb_targets; i++) {
	// 	file >> _risk_thold_ratio[i];
	// }
}


DataHandler::~DataHandler() {
	delete[] _radii;
	delete[] _bdg_rewards_ratio;
	// delete[] _risk_thold_ratio;
	delete[] _target_locs;
}

void DataHandler::print() {
	if (_nb_clusters != 0) {
		cout << "0. targets are clustered: " << _nb_clusters << endl;
	}	
	cout << "1. start depot:  " << _depot1_loc._x << " " << _depot1_loc._y << endl;
	cout << "2. end depot:  " << _depot2_loc._x << " " << _depot2_loc._y << endl;
	cout << "3. # of cirlces:  " << _nb_targets << endl;
	cout << "4. the lengths of radius are:  ";
	for (int i = 0; i < _nb_targets; i++)
		cout << _radii[i] << "  ";

	cout << "\n";
	cout << "5. coordinates of targets' locations are:  \n";
	for (int i = 0; i < _nb_targets; i++)
		cout << _target_locs[i]._x << "\t" << _target_locs[i]._y << "\n";

	cout << "6. lower bounds of rewards:  ";
	for (int i = 0; i < _nb_targets; i++)
		cout << _bdg_rewards_ratio[i] << "  ";
	cout << "\n";
	// cout << "7. upper bounds of risks:  ";
	// for (int i = 0; i < _nb_targets; i++)
	// 	cout << _risk_thold_ratio[i] << "  ";
	// cout << "\n";
}
