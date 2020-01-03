
#ifndef _SPFORMULATION_H_
#define _SPFORMULATION_H_

#include <fstream>
#include <iostream>
#include "gurobi_c++.h"
#include "costInfo.h"

using namespace std;

class SPFormulation{
	
	private:		
		int m, n;
		GRBVar*** x;
		costInfo* cost;
		
	public:
		SPFormulation(costInfo* cost, GRBModel* model);
		void solve(GRBModel* model);
		void printSol(GRBModel* model);
};

#endif


