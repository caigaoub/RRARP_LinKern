#include "SubtourCuts.h"
#include <tuple>
#include <assert.h>

SubtourCuts::SubtourCuts(GRBVar** x_,  int _nb_targets){
	this->_var_x = x_;
	this->_nb_targets = _nb_targets;
	this->_size_var_x = _nb_targets +2;
}


void SubtourCuts::callback() {
	try {
		if (where == GRB_CB_MIPSOL) {
			double **x_sol = new double*[_nb_targets+2];
			for (int i = 0; i < _nb_targets+2; i++) {
				x_sol[i] = new double[_nb_targets+2];
				x_sol[i] = getSolution(_var_x[i], _nb_targets+2);
			}
			// for (int i = 0; i < _size_var_x; i++) {
			// 	for (int j = 0; j < _size_var_x; j++) {
			// 		if(abs(x_sol[i][j]) < 1e-6){
			// 				cout << 0 << "  ";
			// 		}
			// 		if(abs(x_sol[i][j]-1) < 1e-6){
			// 			cout << 1 << "  ";
			// 		}
			// 		// cout << sol[i][j] << "   ";	
			// 	}
			// 	cout << endl;
			// }
			// cout << "-----------------------------------" << endl;
			
			int *tour = new int[_nb_targets+2];
			int len;
			findsubtour(_nb_targets+2, x_sol, &len, tour);
			if (len < _nb_targets+2) {
				GRBLinExpr expr = 0;
				for (int i = 0; i < len; i++) {
					for (int j = i + 1; j < len; j++) {
						expr += _var_x[tour[i]][tour[j]] + _var_x[tour[j]][tour[i]];
					}
				}
				addLazy(expr <= len - 1);
				_CB_nb_subtour_cuts++;
			}

			for (int i = 0; i < _nb_targets+2; i++)
				delete[] x_sol[i];
			delete[] x_sol;
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



void SubtourCuts::findsubtour(int  n, double** sol, int*  tourlenP, int*  tour) {
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
			if (i == n) { // all adj(node) are visited
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

void SubtourCuts::print_xSol(double ** x_sol) {
	for (int i = 0; i < _nb_targets+2; i++) {
		for (int j = 0; j < _nb_targets+2; j++) {
			cout << x_sol[i][j] << "   ";
		}
		cout << endl;
	}

}


SubtourCuts::~SubtourCuts(){
	;
}
