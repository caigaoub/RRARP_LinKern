
#ifndef _MCFORMULATION_H_
#define _MCFORMULATION_H_

#include <fstream>
#include <iostream>
#include "gurobi_c++.h"
#include "costInfo.h"

using namespace std;

class MCFormulation{
	
	private:		
		int m, n;
		GRBVar** x;
        GRBVar*** y;
		costInfo* cost;
		
	public:
		MCFormulation(costInfo* cost, GRBModel* model);
		void solve(GRBModel* model);
		void printSol(GRBModel* model);
};

#endif


