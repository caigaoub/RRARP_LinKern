#include <cmath>
#include <sstream>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <math.h>
#include "gurobi_c++.h"
#include <random>


using namespace std;
#define PI  atan(1.0) * 4

double distance(double *x, double *y, int i, int j)
{
	double dx = x[i] - x[j];
	double dy = y[i] - y[j];
	return sqrt(dx*dx + dy*dy);
}

string itos(int i) { stringstream s; s << i; return s.str(); };


int CreateGraphs(int n, int myMAX)
{

	double depotX = (double)(1 + rand() % myMAX);
	double depotY = (double)(1 + rand() % myMAX);

	double * x = new double[n];
	double * y = new double[n];
	
	double ** weight = new double*[n];
	
	int i;
	for (i = 0; i < n; i++)	{
		weight[i] = new double[n];
	}
	for (i = 0; i < n; i++)	{
		x[i] =(double)( 1 + rand() % myMAX);
		y[i] =(double)( 1 + rand() % myMAX);
	}
	int j;

	for (i = 0; i < n; i++)	{
		for (j = 0; j < n; j++)	{
			if (i != j && x[i] == x[j] && y[i] == y[j])	{
				cout << "ERROR: two points have the same coordinates!" << endl;
				cout << "random again!!" << endl;
				return -1;
			}
		}
	}

	for (i = 0; i < n; i++)	{
		for (j = 0; j < n; j++)	{
			if (i != j)	{
				weight[i][j] = distance(x, y, i, j);
			}
			else {
				weight[i][j] = INFINITY;
			}
		}
	}

	GRBEnv* env = new GRBEnv();
	GRBModel model = GRBModel(*env);

	GRBVar* R = new GRBVar[n];
	for (i = 0; i < n; i++) {
		R[i] = model.addVar(0.0, INFINITY, 1.0, GRB_CONTINUOUS, "R_" + itos(i));
	}
	model.update();

	GRBVar* z = new GRBVar;
	*z = model.addVar(0.0, INFINITY, 1.0, GRB_CONTINUOUS, "z");
	GRBLinExpr expr_obj = *z;

	for (i = 0; i < n; i++) {
		expr_obj +=  1.0 * R[i];
	}

	model.setObjective(expr_obj, GRB_MAXIMIZE);
	model.update();


	// constraints
	GRBLinExpr expr = 0;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			if (i != j) {
				expr = R[i] + R[j] - weight[i][j];
				model.addConstr(expr <= 0, "T1_" + itos(i) + "_" + itos(j));
			}
		}
	}
	model.update();

	for (i = 0; i < n; i++) {
		expr = *z - R[i];
		model.addConstr(expr <= 0, "T2_" + itos(i));
	}
	model.update();

	model.optimize();

	if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL){

		vector<double> radius(n);
		for (i = 0; i < n; i++) {
			radius[i] = R[i].get(GRB_DoubleAttr_X);
		}

		double z_val = z->get(GRB_DoubleAttr_X);

		vector<double> minReward(n, 0.00001);
		vector<double> maxRisk(n, 1000000);

		random_device rd;
		mt19937 gen(rd());
		srand(gen());
		double v1 = (double)rand() / RAND_MAX;
		double v2 = (double)rand() / RAND_MAX;
		if (v1 > v2) {
			double temp = v2;
			v2 = v1;
			v1 = temp;
			
		}
		uniform_real_distribution<> dis(v1, v2);
		vector<double> coefRadius(n);
		for (i = 0; i < n; i++)	{
			coefRadius[i] = dis(gen);
		}

		string s = "sample_n_" + itos(n) + ".txt";
		ofstream myfile(s);

		myfile << n << '\n';
		myfile << depotX << '\t' << depotY << '\n';
		for (i = 0; i < n; i++)
		{
			myfile << radius[i] * coefRadius[i] << '\t';
		}
		myfile << '\n';

		for (i = 0; i < n; i++)
		{
			myfile << minReward[i] << '\t';
		}
		myfile << '\n';

		for (i = 0; i < n; i++)
		{
			myfile << maxRisk[i] << '\t';
		}
		myfile << '\n';

		for (i = 0; i < n; i++)
		{
			myfile << x[i] << '\t' << y[i] << '\n';
		}

		myfile.close();
	}


	delete x;
	delete y;
	for (i = 0; i < n; i++)
	{
		delete [] weight[i];
	}
	delete [] weight;

	delete z;
	delete[] R;

}