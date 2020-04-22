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
		_x = other._x;
		_y = other._y;
		return *this;
	}
	inline void print(){cout << "Vertex: (" << _x << "," << _y << ")" << endl;};
};

struct ObsPoint {
	double _rad;
	double _angle;
	double _rewVal;
	ObsPoint& operator= (ObsPoint other) {
		_rad = other._rad;
		_angle = other._angle;
		_rewVal = other._rewVal;
		return *this;
	}
	inline void print(){cout << "Observation point: (" << _rad << "," << _angle << "," << _rewVal << ")" << endl;};
};

class DataHandler {
public:
	string          _name;
	int				_nb_targets;
	int 			_nb_obps;	
	Vertex			_depot1_loc;
	Vertex			_depot2_loc;
	Vertex*			_target_locs;
	double*			_radii;
	// double*			_bdg_rewards_ratio; // budgetary rewards ratio
	vector<vector<ObsPoint>> _all_obps;
	vector<vector<vector<double>>> _all_innermatr;

	DataHandler() {};
	~DataHandler();

	void parse(string filename);
	void read_instance(string);
	void pick_trigraphs(string);

	

	inline double stod(string str){double d; stringstream iss(str); return iss >>d;};
	void print();
};



#endif // !_DATAHANDLER_H_
