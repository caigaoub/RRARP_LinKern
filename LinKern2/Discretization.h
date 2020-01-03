#ifndef _DISCRETIZATION_H_
#define _DISCRETIZATION_H_


#include <vector>
#include <iostream>
#include <string>
#include "DataHandler.h"
#include "CostInfo.h"

using namespace std;

// Coordinates of points on boundaries
struct coord
{
	double x;
	double y;
};

// --
struct VERTEX {
	int tarID;
	int subKey;
	int matKey;
};

struct INDEX {
	int tarID; // target
	int subKey; // discritized point
};

class Discretization
{
private:
	int num_dstzn;
	int num_targets;
	double para_ext_risk;
	double para_int_risk;

	// receive parameters
	double depot_x;
	double depot_y;
	double * center_x;
	double * center_y;
	vector<double>  radii;
	vector<double> min_reward_bounds;
	vector<double> max_risk_bounds;

	CostInfo* CI;
	double subarc_angle;
	int	num_nodes_sglIdx;

	vector<vector<coord>> nodes;

	vector<vector<double>> intDist;
	vector<vector<double>> sglIdxMat;

	vector<int> num_subarc_chord;

public:
	Discretization(int k, DataHandler *, CostInfo*);
	~Discretization();

	inline int get_num_circles() { return num_targets; }
	inline int get_num_dstzn() { return num_dstzn; }
	inline double get_depot_x() { return depot_x; }
	inline double get_depot_y() { return depot_y; }

	inline vector<vector<double>>& getExtDist() { return sglIdxMat; };
	inline vector<vector<double>>& getIntDist() { return intDist; };

	void generate_nodes_crds();
	void cal_ext_dist();
	void cal_int_dist();
//	void generate_sglIdxMat();

	double evaluate_risk(double radius, double w);
	double evalute_reward(double radius, double w);
	
	void build();
	
//	void fast_sp(vector<vector<double>> & SDS, vector<int> &seq);
	void fast_sp2(vector<VERTEX>& Tour, double& optimal, vector<int> &seq);

};

#endif // !_DISCRETIZATION_H_
