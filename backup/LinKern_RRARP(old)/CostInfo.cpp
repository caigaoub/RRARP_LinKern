#include <iostream>
#include <math.h>
#include "CostInfo.h"
#include "DataHandler.h"
#include <cmath>

// CostInfo
CostInfo::CostInfo(DataHandler *graph) {
	
	// N is the number of nodes in TSP
	this->N = graph->get_num_targets() + 2;

	X = new double[N];
	Y = new double[N];

	min_dist_nbds = new double*[N];

	X[0] = graph->get_depot_x();
	Y[0] = graph->get_depot_y();

	X[N - 1] = graph->get_depot_x();
	Y[N - 1] = graph->get_depot_y();
	int i;
	for (i = 1; i <= N - 2; i++)
	{
		X[i] = graph->get_center_x()[i - 1];
		Y[i] = graph->get_center_y()[i - 1];
	}

	for (i = 0; i < N; i++)
	{
		min_dist_nbds[i] = new double[N];
	}

	for ( i = 1; i < N - 1; i++)
	{
			min_dist_nbds[0][i] = sqrt(pow(X[0] - X[i], 2) + pow(Y[0] - Y[i], 2)) - graph->get_radii()[i - 1];
			min_dist_nbds[i][0] = min_dist_nbds[0][i];
			min_dist_nbds[N - 1][i] = min_dist_nbds[0][i];
			min_dist_nbds[i][N - 1] = min_dist_nbds[0][i];
	}
	min_dist_nbds[0][N - 1] = 0.0;
	min_dist_nbds[N - 1][0] = 0.0;


	for (int i = 1; i < N - 1; i++)
	{
		for (int j = i + 1; j < N - 1; j++)
		{	
			min_dist_nbds[i][j] = sqrt(pow(X[i] - X[j], 2) + pow(Y[i] - Y[j], 2)) - graph->get_radii()[i - 1] - graph->get_radii()[j - 1];
			min_dist_nbds[j][i] = min_dist_nbds[i][j];
		}
	}
	//diagonal
	for (int i = 0; i < N; i++)
	{
		min_dist_nbds[i][i] = 0.0;
	}
	/*
	
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			cout << min_dist_nbds[i][j] << "  ";
		}
		cout << "\n";
	}
	*/
	
}

void CostInfo::printCostInfo()
{
	cout << "--------------- The weights of edges for TSP---------------" << endl;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			cout << min_dist_nbds[i][j] << '\t' << " ";
		}
		cout << endl;
	}
	cout << "-------------- Fancy Cutting Line -------------------" << endl;
	cout << endl;
}

CostInfo::~CostInfo() {
	for (int i = 0; i < N; i++)
	{
		delete(min_dist_nbds[i]);
	}
	delete X;
	delete Y;

}


