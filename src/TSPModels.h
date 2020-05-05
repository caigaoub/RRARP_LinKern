#ifndef _TSPMODEL_H_
#define _TSPMODEL_H_
#include "gurobi_c++.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>
// #include "BendersCuts.h"
using namespace std;

class TSPModel_STE // subtour elimination formulation
{
public:
	double **                   _cost;
	bool 						_edge_cost_status = false;

	GRBModel* 					_model;
	int                         _size_var_x;
	GRBVar**                    _var_x;
	int                         _nb_subtour_cuts = 0;
	vector<int>					_opt_seq;
	double 						_curbestobj = -911;

	TSPModel_STE(){};
	~TSPModel_STE();
	void init_edge_weights(int, const vector<vector<double>> &);
	void create_formula();
	void solve();
//	void writeSol(GRBModel * model);
	// void printSol(GRBModel* model);
	inline static string itos(int i) {stringstream s; s << i; return s.str();};
    inline static string dtos(double i) {stringstream s; s << i; return s.str();};
};

	


class SubtourElimCuts : public GRBCallback
{
public:
	int 				_size_var_x;
	GRBVar ** 			_var_x;
	int 				_nb_subtour_cuts = 0;

	SubtourElimCuts(GRBVar** varX_arg, int _size_varX_arg){
		this->_var_x = varX_arg;
		this->_size_var_x = _size_varX_arg;
	};
	static void find_subtour(int  n, double** sol, int*  tourlenP, int*  tour){
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
	};

protected:
	void callback(){
		try {
			if (where == GRB_CB_MIPSOL)	{
				double **x = new double*[_size_var_x];
				int *tour = new int[_size_var_x];
				int len;
				for (int i = 0; i < _size_var_x; i++){
					x[i] = new double[_size_var_x];
					x[i] = getSolution(_var_x[i], _size_var_x);
				}
				find_subtour(_size_var_x, x, &len, tour);

				if (len < _size_var_x){
					// Add subtour elimination constraint
					GRBLinExpr expr = 0;
					for (int i = 0; i < len; i++)
						for (int j = i + 1; j < len; j++)
							expr += _var_x[tour[i]][tour[j]] + _var_x[tour[j]][tour[i]];
					addLazy(expr <= len - 1);
					_nb_subtour_cuts++;
				}
				for (int i = 0; i < _size_var_x; i++)
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
	};
};


#endif // !_TSPMODEL_H_
