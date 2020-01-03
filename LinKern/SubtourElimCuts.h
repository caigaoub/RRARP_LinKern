#ifndef _SUBTOURELIM_H_
#define _SUBTOURELIM_H_

#include "gurobi_c++.h"
#include <cstdlib>
#include <cmath>
#include <sstream>
#include "costInfo.h"
#include <algorithm>
using namespace std;


class SubtourElimCuts : public GRBCallback
{
private:

	int N;
	GRBVar ** y;

	unsigned int numSubtourCuts;

public:
	SubtourElimCuts(GRBVar** yVars, int yn);
	inline int get_numSubtourCuts() { return numSubtourCuts; }

protected:
	void callback();

};

#endif // !_SUBTOURELIM_H_
