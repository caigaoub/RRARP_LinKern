#ifndef _STEFORMULATION_H_
#define  _STEFORMULATION_H_

#include <sstream>
#include "gurobi_c++.h"
#include "DataHandler.h"
#include "SubtourCuts.h"
// #include "ProgTime.h"
#include <boost/filesystem.hpp>
#include <sys/stat.h>
#include "MacroName.h"

using namespace std;

class STEFormulation {
public:
	DataHandler *												_dataset = nullptr;

	int 														_graphsize;
	vector<vector<pair<bool,double>>> * 						_G_dic;
	vector<vector<vector<pair<double, string>>>>* 				_all_innerG;
	int 														_nb_dstzn;
	int 														_nb_targets;


	GRBModel*													_model;
	int															_size_var_x;
	GRBVar**													_var_x;
	int															_size_var_u;
	GRBVar**													_var_u;

	int															_total_nb_subtour_cuts = 0;



	STEFormulation() {};
	~STEFormulation();
	void build_formul(GRBModel*, DataHandler&);
	
	void solve_formul_wCB();// with callback

	void get_optimal_sol(double **);
	string itos(int i) { stringstream s; s << i; return s.str(); }
	// void write_solution(string, int);
};





#endif
