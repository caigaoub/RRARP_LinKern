#pragma once

#include <vector>
#include "Discretization.h"
#include "CostInfo.h"
using namespace std;

struct LKEdgeX {
	double dist; // distance
	int spKeyA; 
	int spKeyB;
	int spKeyC;
	int spKeyD;	
};

struct  LKEdgeY {
	double dist;
	VERTEX newPointA; // circle with smaller key
	VERTEX newPointB;
};

struct TARGET {
	int id;
	int key;
};

//typedef pair<double, pair<int, int>> TwoEndpoints; // distance, pi, pi+1

class LinKern {
private: 
	Discretization* disc;

	int N;
	int num_dstzn;

	vector<vector<coord>> nodes;

	vector<vector<double>>  D;
	vector<vector<double>> intDist;
	

	vector<VERTEX> optimalSPTour;
	double optimalLens;
//	vector<vector<pair<double, INDEX>>> sortedCost;
	int runTimes;

public:
	LinKern(Discretization*, CostInfo*);
	~LinKern();

	void solve();
	void kingsLanding();

	double LKSearch(vector<VERTEX>&, double&);
	double LKSearch(vector<VERTEX> &, int);
	double LKSearch(vector<VERTEX>&, int, int);
	double LKSearch_x2(vector<VERTEX>&, vector<VERTEX>&, vector<pair<int, int>>&, vector<vector<VERTEX>>&, int&, double&, double&, pair<VERTEX, VERTEX>&, double);
	double LKSearch_Continue(vector<VERTEX>&, vector<VERTEX>&, vector<pair<int, int>>&, vector<vector<VERTEX>>&, int&, double&, double&, pair<VERTEX, VERTEX>&, double);


	bool isPositiveGain(double Gi, double g_local);

	bool updateGStar(double&, double, double, double);

	int getClosingUpNode(vector<VERTEX>&, vector<pair<TARGET, TARGET>>&, int);

	LKEdgeX get_x(vector<VERTEX>&, int, int);
	LKEdgeY get_y(vector<VERTEX>&, VERTEX, VERTEX);
	LKEdgeY get_y_type1(vector<VERTEX>&, VERTEX, VERTEX);
	LKEdgeY get_y_type2(vector<VERTEX>&, VERTEX, VERTEX);
	LKEdgeY get_y_type3(vector<VERTEX>&, VERTEX, VERTEX);

	bool isInSet(vector<pair<int, int>>&, int);
	
	void updatePoints(vector<VERTEX>&, LKEdgeX, LKEdgeX, LKEdgeY, LKEdgeY);
	
//	void FLIP(vector<VERTEX>&, vector<pair<TARGET, TARGET>>, int k);
	void flip(vector<VERTEX>&, int a, int b, int c, int d);


	vector<int> createRandomTour();
	void swap(vector<int>&, int, int);
	
	int getSPKey(vector<VERTEX>&, VERTEX);
	int getKey(vector<VERTEX>&, int);
	int trueKey(int);
	int matrixKey(int, int, bool);
	bool getDirection(int, int);
	bool isNeighbor(vector<VERTEX>&, int, int);

	double getTourLens(vector<VERTEX>&);
	vector<int> getTspTour(vector<VERTEX>&);
};
