#ifndef _DATAHANDLER_H_
#define _DATAHANDLER_H_

#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <limits>
using namespace std;

struct Vertex {
	double _x;
	double _y;
	Vertex& operator= (Vertex other) {
		// swap(_x, other._x);
		// swap(_y, other._y);
		_x = other._x;
		_y = other._y;
		return *this;
	}
	inline void print(){cout << "Vertex: (" << _x << "," << _y << ")" << endl;};
};

class DataHandler {
public:
	string          _name;
	int				_nb_targets;
	int				_nb_clusters = 0; // default no cluster
	Vertex			_depot1_loc;
	Vertex			_depot2_loc;
	Vertex*			_target_locs;
	double*			_radii;
	double*			_bdg_rewards_ratio; // budgetary rewards ratio
	// double*			_risk_thold_ratio; // maximum risk taken at that inner trajectory	

	DataHandler() {};
	~DataHandler();
	void parse(string filename);
	void print();
};

#endif // !_DATAHANDLER_H_
