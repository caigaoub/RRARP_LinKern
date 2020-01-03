#ifndef _STEFORMULATION_H_
#define _STEFORMULATION_H_


#include "gurobi_c++.h"
#include "costInfo.h"
#include "SubtourElimCuts.h"
using namespace std;

class STEFormulation
{
private:
	int N;

	GRBVar** y;

	costInfo * cost;


	int numSubtourConst;
	unsigned int status;
public:
	STEFormulation(costInfo * c, GRBModel * model);
	void solve(GRBModel * model);
	void printSol(GRBModel* model);


	string itos(int i) { stringstream s; s << i; return s.str();}
};
#endif // !_STEFORMULATION_H_
