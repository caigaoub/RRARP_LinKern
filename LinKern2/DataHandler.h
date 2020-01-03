#ifndef _DATAHANDLER_H_
#define _DATAHANDLER_H_

#include <iostream>
#include <string>
#include <fstream>
#include <set>
#include <vector>

using namespace std;
struct Point {
	double x;
	double y;
};

struct LineSegment{
	Point endpoint_1;
	Point endpoint_2;
};

struct TangentLineSegment {
	LineSegment OuterTL1;
	LineSegment OuterTL2;
	LineSegment InnerTL1;
	LineSegment InnerTL2;
};



class DataHandler
{
private:
	int			 num_targets;
	double		 depot_x;
	double		 depot_y;
	double*		 radii;			// length of radius for each circle
	double*		 min_reward_bounds;	
	double*		 max_risk_boounds;
	double*		 center_x;      // x-coordinate of the center
	double*		 center_y;      // y-coordinate of the center

public:
	DataHandler(const char* filename);

	~DataHandler();

	inline int      get_num_targets() { return num_targets; }
	inline double   get_depot_x() { return depot_x; }
	inline double   get_depot_y() { return depot_y; }
	inline double*  get_center_x() { return center_x; }
	inline double*  get_center_y() { return center_y; }
	inline double*  get_radii() { return radii; }
	inline double*  get_min_rewards() { return min_reward_bounds; }
	inline double*  get_max_risks() { return max_risk_boounds; }

	inline void swap(double& a, double& b) { double c = b; b = a; a = c; }
	void getTangentLines(int, int);
	void printGraphInfo();

};

#endif // !_DATAHANDLER_H_

