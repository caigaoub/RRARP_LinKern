
#ifndef _MTZFORMULATION_H_
#define _MTZFORMULATION_H_

#include <fstream>
#include <iostream>
#include "gurobi_c++.h"
#include "costInfo.h"

using namespace std;

class MTZFormulation{
	
	private:		
		int m, n;
		GRBVar** x;
        GRBVar* u;
		costInfo* cost;
		
	public:
		MTZFormulation(costInfo* cost, GRBModel* model);
		void solve(GRBModel* model);
		void printSol(GRBModel* model);
};

#endif


