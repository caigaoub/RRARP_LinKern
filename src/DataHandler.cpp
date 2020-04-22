#include "DataHandler.h"



void DataHandler::read_instance(string filename) {
	fstream file(filename);
	if (!file) {
		cerr << "ERROR: could not open file '" << filename << "' for reading'" << endl;
		throw(-1);
	}
	// Geographic Information of Targets
	auto pos = filename.find_last_of("/");
    _name = filename.substr(pos+1, filename.size());
    // size_t found = _name.find("c");
    double temp;
	file >> _nb_targets >> temp;
	file >> _depot1_loc._x >> _depot1_loc._y;
	file >> _depot2_loc._x >> _depot2_loc._y;
	_target_locs = new Vertex[_nb_targets];
	_radii = new double[_nb_targets];
	for (int i = 0; i < _nb_targets; i++) {
		file >> _target_locs[i]._x >> _target_locs[i]._y >> _radii[i];
	}

	// _bdg_rewards_ratio = new double[_nb_targets];
	// _risk_thold_ratio = new double[_nb_targets];
	// for (int i = 0; i < _nb_targets; i++) {
	// 	file >> _bdg_rewards_ratio[i];
	// }
}


void DataHandler::pick_trigraphs(string dir){

	vector<int> pickarr(_nb_targets, 0);
	for(unsigned i=0; i < _nb_targets; i++)
		pickarr[i] = i+1;

	_all_innermatr.resize(_nb_targets);
	for(int t=0; t<_nb_targets; t++){
		_all_innermatr[t].resize(16);
		for(int i=0; i<16; i++){
			_all_innermatr[t][i].resize(16, 0);
		}		
	}
	for(auto itr = pickarr.begin(); itr != pickarr.end(); itr++){
		string filename = dir + to_string((*itr)) + ".trigraph_k16";
		fstream file(filename);
		if (!file) {
			cerr << "ERROR: could not open file '" << filename << "' for reading'" << endl;
			throw(-1);
		}
		// file >> _nb_obps >> 



	}
	_all_obps.resize(_nb_targets);
	// _nb_obps  = 100;
	vector<string> emptystr(_nb_obps);
	string::size_type sz = 0;
	ObsPoint obp; 
	for (int i = 0; i < _nb_targets; i++){
		for(int j = 0; j < _nb_obps; j++){
			file >> emptystr[j];
		}
		_all_obps[i].resize(_nb_obps);
		for(int j = 0; j < _nb_obps; j++){
			pos = emptystr[j].find_first_of(":");
			string temp = emptystr[j].substr(0, pos);
			obp._rad = stod(temp.substr(sz));
			string temp_left = emptystr[j].substr(pos+1, emptystr[j].size());
			pos = temp_left.find_first_of(":");
			temp = temp_left.substr(0, pos);
			obp._angle = stod(temp.substr(sz));
			temp = temp_left.substr(pos+1, temp_left.size());
			obp._rewVal = stod(temp.substr(sz));
			// obp.print();
			_all_obps[i][j] = obp;
		}	
	}
}


void DataHandler::parse(string filename) {
	fstream file(filename);
	if (!file) {
		cerr << "ERROR: could not open file '" << filename << "' for reading'" << endl;
		throw(-1);
	}
	// Geographic Information of Targets
	auto pos = filename.find_last_of("/");
    _name = filename.substr(pos+1, filename.size());
    // size_t found = _name.find("c");
	// if (found == string::npos) {
	file >> _nb_targets >> _nb_obps;
	_nb_obps += 16;
	// }
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

	_all_obps.resize(_nb_targets);
	vector<string> emptystr(_nb_obps);
	string::size_type sz = 0;
	ObsPoint obp; 
	for (int i = 0; i < _nb_targets; i++){
		for(int j = 0; j < _nb_obps; j++){
			file >> emptystr[j];
		}
		_all_obps[i].resize(_nb_obps);
		for(int j = 0; j < _nb_obps; j++){
			pos = emptystr[j].find_first_of(":");
			string temp = emptystr[j].substr(0, pos);
			obp._rad = stod(temp.substr(sz));
			string temp_left = emptystr[j].substr(pos+1, emptystr[j].size());
			pos = temp_left.find_first_of(":");
			temp = temp_left.substr(0, pos);
			obp._angle = stod(temp.substr(sz));
			temp = temp_left.substr(pos+1, temp_left.size());
			obp._rewVal = stod(temp.substr(sz));
			// obp.print();
			_all_obps[i][j] = obp;
		}	
	}
}


DataHandler::~DataHandler() {
	delete[] _radii;
	delete[] _bdg_rewards_ratio;
	delete[] _target_locs;
}

void DataHandler::print() {
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
