#include "SubtourElimCuts.h"

SubtourElimCuts::SubtourElimCuts(GRBVar** yVars, int yn)
{
	y = yVars;
	N = yn;

	numSubtourCuts = 0;

}


void findsubtour(int  n, double** sol, int*  tourlenP, int*  tour)
{
	bool* seen = new bool[n];
	int bestind, bestlen;
	int i, node, len, start;

	for (i = 0; i < n; i++)
		seen[i] = false;

	start = 0;
	bestlen = n + 1;
	bestind = -1;
	node = 0;
	while (start < n) {
		for (node = 0; node < n; node++)
			if (!seen[node])
				break;
		if (node == n)
			break;
		for (len = 0; len < n; len++) {
			tour[start + len] = node;
			seen[node] = true;
			for (i = 0; i < n; i++) {
				if (sol[node][i] > 0.5 && !seen[i]) {
					node = i;
					break;
				}
			}
			if (i == n) {
				len++;
				if (len < bestlen) {
					bestlen = len;
					bestind = start;
				}
				start += len;
				break;
			}
		}
	}

	for (i = 0; i < bestlen; i++)
		tour[i] = tour[bestind + i];
	*tourlenP = bestlen;

	delete[] seen;
}

void SubtourElimCuts::callback()
{
	try {
		if (where == GRB_CB_MIPSOL) 
		{
			// Found an integer feasible solution - does it visit every node?
			double **x = new double*[N];
			int *tour = new int[N];
			int i, j, len;
			for (i = 0; i < N; i++)
			{
				x[i] = new double[N];
				x[i] = getSolution(y[i], N);
			}
		//	cout << "before: " << endl;
		//	for (i = 0; i < N; i++)
			{
			//	for (j = 0; j < N; j++)
			//	{
					
				//	cout << x[i][j] << '\t';
			//	x[i] = getSolution(y[i], N);
				//}
			//	cout << '\n';
			}
			/*
			for (i = 0; i < N; i++)
			{
				for (j = 0; j < N; j++)
				{
					x[i][j] = getSolution(y[i][j]);
					if (x[i][j] > 0.5)
					{
						x[j][i] = x[i][j];
					}
				}
			}
			cout << "after: " << endl;
			for (i = 0; i < N; i++)
			{
				for (j = 0; j < N; j++)
				{
					cout << x[i][j] << '\t';
				}
				cout << '\n';
			}
			*/
			findsubtour(N, x, &len, tour);

			if (len < N) 
			{
				// Add subtour elimination constraint
				GRBLinExpr expr = 0;
				for (i = 0; i < len; i++)
					for (j = i + 1; j < len; j++)
						expr += y[tour[i]][tour[j]] + y[tour[j]][tour[i]];
				addLazy(expr <= len - 1);

				numSubtourCuts++;
			}
			for (i = 0; i < N; i++)
				delete[] x[i];
			delete[] x;
			delete[] tour;
		}
	}
	catch (GRBException e) {
		cout << "Error number: " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during callback" << endl;
	}
}

