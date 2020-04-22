#ifndef _LKALGO_H_
#define _LKALGO_H_


#include "DataHandler.h"
#include <tuple>
#include <cmath>
#include <math.h>
#include <deque>
#include <vector>
#include <numeric>
#include "ProgTime.h"
#include <string>
#define INF numeric_limits<double>::infinity()


using namespace std;

class


class LK{
public:
	int 									_graphsize;
	vector<vector<double>>					_riskgraph;
	vector<vector<vector<double>>>			_innerriskgraph;

	// vector<int> 

	void gain_search();

};













#endif 