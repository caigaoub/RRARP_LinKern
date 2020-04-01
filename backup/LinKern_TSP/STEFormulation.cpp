/* Copyright 2016, Gurobi Optimization, Inc. */

/* Solve a traveling salesman problem on a randomly generated set of
points using lazy constraints.   The base MIP model only includes
'degree-2' constraints, requiring each node to have exactly
two incident edges.  Solutions to this model may contain subtours -
tours that don't visit every node.  The lazy constraint callback
adds new constraints to cut them off. */

#include "gurobi_c++.h"
#include <cstdlib>
#include <cmath>
#include <sstream>
#include "STEFormulation.h"
#include "costInfo.h"

using namespace std;

extern void findsubtour(int  n, double** sol, int*  tourlenP, int*  tour);


STEFormulation::STEFormulation(costInfo* c, GRBModel * model)
{
	// initialize some data
	
	this->cost = c;
	N = cost->getNodes();
	
	// set up the model
	model->set(GRB_IntAttr_ModelSense, 1);

	// Must set LazyConstraints parameter when using lazy constraints
	model->getEnv().set(GRB_IntParam_LazyConstraints, 1);
	model->getEnv().set(GRB_IntParam_PreCrush, 1);
	model->getEnv().set(GRB_IntParam_MIPFocus, 1); // focus on feasibility

	// add variables
	y = new GRBVar*[N];

	int i, j;
	for (i = 0; i < N; i++)
		y[i] = new GRBVar[N];

	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			y[i][j] = model->addVar(0.0, 1.0, cost->getCost(i,j), GRB_BINARY, "y_" + itos(i) + "_" + itos(j));

			//y[j][i] = y[i][j];
		}
	}
	model->update();

	// Degree-1 constraints
	for (i = 0; i < N; i++)
	{
		GRBLinExpr expr = 0;
		for (j = 0; j < N; j++)
		{
			expr += 1.0* y[i][j];
		}
		model->addConstr(expr == 1, "deg1_row" + itos(i));
		GRBLinExpr expr2 = 0;
		for (j = 0; j < N; j++)
		{
			expr2 += 1.0* y[j][i];
		}
		model->addConstr(expr2 == 1, "deg1_col" + itos(i));
	}

	// Forbid edge from node back to itself
	for (i = 0; i < N; i++)
	{
		y[i][i].set(GRB_DoubleAttr_UB, 0);
	}

	model->update();
}

void STEFormulation::solve(GRBModel *model)
{

	try
	{

		SubtourElimCuts cb = SubtourElimCuts(y, N);
		model->setCallback(&cb);
		// Optimize model
		model->optimize();
		numSubtourConst = cb.get_numSubtourCuts();
		cout << numSubtourConst << endl;
		if (model->get(GRB_IntAttr_Status) == GRB_OPTIMAL)
			status = 0;

	}
	catch (GRBException e) {
		cout << "Error number: " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;

	}
	catch (...) {
		cout << "Error during optimization" << endl;
	}

}

void STEFormulation::printSol(GRBModel *model)
// Extract solution
{
	if (model->get(GRB_IntAttr_SolCount) > 0) 
	{
		double **sol = new double*[N];
		int i;
		for (i = 0; i < N; i++)
			sol[i] = model->get(GRB_DoubleAttr_X, y[i], N);

		int* tour = new int[N];
		int len;

		findsubtour(N, sol, &len, tour);

		cout << "Tour: ";
		for (i = 0; i < len; i++)
			cout << tour[i] << " ";
		cout << endl;

		for (i = 0; i < N; i++)
			delete[] sol[i];
		delete[] sol;
		delete[] tour;

	}
	for (int i = 0; i < N; i++)
		delete[] y[i];
	delete[] y;
	
}