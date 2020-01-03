#include <vector>
#include "gurobi_c++.h"
#include "DataHandler.h"
#include "Discretization.h"
#include "CostInfo.h"
#include "LinKern.h"
#include <ctime>
using namespace std;
extern int CreateGraphs(int n, int myMAX);

void main22222()
{
	CreateGraphs(15, 220);
}
void main()
{
	try
	{
		DataHandler graph("sample_n_13_v_5.txt");
		CostInfo cost(&graph);

		int k = 8;
		Discretization disc(k, &graph, &cost);

		clock_t begin = clock();

		LinKern lkObj(&disc, &cost);
		lkObj.solve();
		
		clock_t end = clock();
		double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

		cout << "Time Elapsed:  " << elapsed_secs << endl;
		
	}

	catch (GRBException& ex)
	{
		cout << "Error number: " << ex.getErrorCode() << endl;
		cout << ex.getMessage() << endl;
	}
	catch (...)
	{
		cerr << "Error" << endl;
	}

	system("pause");
}
