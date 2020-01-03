#ifndef _COSTINFO_H_
#define _COSTINFO_H_

#include <fstream>
#include <sstream>
#include <iostream>
#include "DataHandler.h"


using namespace std;

class CostInfo {

private:
	int     N; // capital N: number of nodes in TSP;
	double*  X; // capital X;
	double*  Y; // capital Y;
	double** min_dist_nbds;

public:
	CostInfo(DataHandler *graph);
	~CostInfo();


	inline int get_num_targets() { return N - 2; };
	inline double ** get_min_dist_nbds() { return min_dist_nbds; }

	void printCostInfo();

	inline static string itos(int i) { stringstream s; s << i; return s.str(); };
	inline static string dtos(double i) { stringstream s; s << i; return s.str(); };
};

#endif