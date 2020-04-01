#pragma once
#include <vector>

#include "costInfo.h"

using namespace std;


struct Vertex {
	int id;
	int key;
};

struct Edge {
	int u;
	int v;
};



class LinKern {
private:
	int N;

	double** D; // distance matrix
	vector<vector<pair<double, int>>> sortedCost;

	vector<int> optimalTour;

	double optimalLength;

	// control parameters
	int runTimes;

public:
	LinKern(costInfo*);
	
	~LinKern();

	void solve();
	void kingsLanding();
	
	double LKSearch(vector<int>& Tour); // return minimal length for current Tour
	pair<double,double> LKSearch(vector<int>&, int);
	pair<double,double> LKSearch(vector<int>&, int, int);

	double LKSearch_x2(vector<int>&, vector<int>&, vector<pair<Vertex, Vertex>>&, int&, double&, double&);
	double LKSearch_Continue(vector<int>&, vector<int>&, vector<pair<Vertex, Vertex>>&, int&, double&, double&);

	double LKSearch_alternative_x2(vector<int>&, vector<pair<Vertex, Vertex>>&, double&, double&);
	pair<vector<int>, vector<int>> getTwoSegments(vector<int>&, int, int, int, int);
	bool isInSet(vector<int>&, int);
//	bool updateGStar(double&, double&,);

	bool isInSet(vector<pair<Vertex, Vertex>>& X, int);
	bool isInTour(vector<int>&, Edge);
	bool isSameEdge(Edge, Edge);
	bool isNeighbor(vector<int>&, int, int);
	bool isPositiveGain(double, double);

	void updateGain(double&, double);
	bool updateGStar(double&, double Gi,  int, int, int);


	int getClosingUpNode(vector<int>&, vector<pair<Vertex, Vertex>>& , int);

	void FLIP(vector<int>&, vector<pair<Vertex, Vertex>> X,  int k);
	void flip(vector<int>&, int a, int b, int c, int d);

	vector<int> createRandomTour();
	void swap(vector<int>&, int, int);
	

	
	double getTourLength(vector<int>);
	void printSol(vector<int>&);
	int trueKey(int);
	int getKey(vector<int>&, int);
};

