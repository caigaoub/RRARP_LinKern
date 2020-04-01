
#ifndef _QAPFORMULATION_H_
#define _QAPFORMULATION_H_

#include <fstream>
#include <iostream>
#include "gurobi_c++.h"
#include "costInfo.h"

using namespace std;

class QAPFormulation{
	
	private:		
		int m, n;
		GRBVar** x;
        GRBVar*** u;
		costInfo* cost;
		
	public:
		QAPFormulation(costInfo* cost, GRBModel* model);
		void solve(GRBModel* model);
		void printSol(GRBModel* model);
};

#endif


