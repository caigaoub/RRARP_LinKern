#ifndef _DATAHANDLER_H_
#define _DATAHANDLER_H_

#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <limits>
#include <tuple>

using namespace std;


struct Vertex {
	double		_x;
	double		_y;
	Vertex& operator= (Vertex other) {
		_x = other._x;
		_y = other._y;
		return *this;
	}
	inline void print(){cout << "Vertex: (" << _x << "," << _y << ")" << endl;};
};

class myVector {
public:
	double _x;
	double _y;
	myVector(Vertex initialp, Vertex endp) {
		this->_x = endp._x - initialp._x;
		this->_y = endp._y - initialp._y;
	}
	~myVector() {};
	inline double get_vecLen() {return sqrt(_x*_x + _y*_y); };
};


struct ObsPoint {
	double		_rad;
	double		_angle;
	double		_rewVal;
	ObsPoint& operator= (ObsPoint other) {
		_rad = other._rad;
		_angle = other._angle;
		_rewVal = other._rewVal;
		return *this;
	}
	inline void print(){cout << "Observation point: (" << _rad << "," << _angle << "," << _rewVal << ")" << endl;};
};

class DataHandler{
public:
	string  												_name;
	int														_nb_targets;
	int 													_nb_obps;	
	Vertex													_depot1_loc;
	Vertex													_depot2_loc;
	vector<Vertex>											_target_locs;
	vector<double>											_radii;
	vector<double>											_bdg_rewards_ratio; // budgetary rewards ratio
	vector<vector<ObsPoint>> 								_all_obps;


	/* risk graph info*/
	vector<vector<Vertex>>									_points; // all turning points after boundary partitioning
	int 													_nb_dstzn = 8;
	int 													_graphsize;
	double 													_par_outer = 1.0;
	vector<vector<pair<bool,double>>>						_riskgraph;
	vector<vector<pair<bool,double>>>						_riskgraph_dic;
	vector<vector<vector<pair<double, string>>>> 			_all_innermatr;
	double 													_total_inrisk_LB = 0;

	vector<vector<double>>									_c2c_adjmatr; // center to center
	vector<vector<double>>									_b2b_adjmatr; // boudary to boundary


	DataHandler() {};
	~DataHandler(){};

	void read_instance(string);
	void pick_trigraphs(string);
	void build_riskgraph();
	double eucl_distance(Vertex&, Vertex&);
	tuple<bool, double,double> is_intersected(Vertex v, Vertex u, int tar);
	double get_risk_outerTrajc(Vertex v, Vertex u);

};



#endif // !_DATAHANDLER_H_
