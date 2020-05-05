#ifndef _SUBTOURCUTS_H_
#define _SUBTOURCUTS_H_

#include <cmath>
#include <stdlib.h>
#include <ctime>
#include <chrono>
#include <queue>
#include <sstream>
#include "gurobi_c++.h"

using namespace std;

class SubtourCuts : public GRBCallback {
private:
  	int 										_nb_dstzn;
	int 										_nb_targets;
	int 										_size_var_x;
	GRBVar ** 									_var_x;
	int 										_CB_nb_subtour_cuts=0;

public:
	~SubtourCuts();
	SubtourCuts(GRBVar**, int);
	static void findsubtour(int, double**, int*, int* );
	inline int get_nb_subtour_cuts() { return _CB_nb_subtour_cuts; }
	inline string itos(int i) { stringstream s; s << i; return s.str(); }
	void print_xSol(double ** x_sol);

protected:
	void callback();
};


#endif // !_SUBTOURCUTS_H_
