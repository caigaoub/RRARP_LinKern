#ifndef _PULSEALGO_H_
#define _PULSEALGO_H_
#include "DataHandler.h"
#include <tuple>
#include <cmath>
#include <math.h>
#include <deque>
#include <vector>
#include <numeric>
#define INF numeric_limits<double>::infinity()


class Pulse{
public:
	DataHandler *					_dataset = nullptr;
	int                             _taridx = -1;
	int 							_graphsize;
	vector<vector<double>> 			_graph;
	vector<vector<bool>>			_tmp_graph;

	double 							_total_reward = 0;
	double 							_budget_pct = -1.0;
	double 							_budget = -1.0;
	int								_source = 0;
	int 							_sink = -1;
	double 							_curbest_objval = INF;
	vector<int>						_curbest_path;
	vector<double> 					_lbrisk_sink;
	vector<vector<int>>				_lbriskpaths_sink;
	vector<double>					_rewards_lbrpaths_sink;		

	vector<double>					_rewards_etd;

	int 							_nb_tnps = -1;
	vector<Vertex>					_tnps_tar;
	int 							_iterations = 0;

	vector<pair<double,double>>		_domi_labels;






	Pulse() {};
	~Pulse() {};

	void initialize_generalgraph(vector<vector<double>>&, vector<double>&, int, double);
	void initialize(DataHandler&, int, int, int);
	void initialize(DataHandler&, int, int, int, string, int, int);

	void update_graph(Vertex&, Vertex&);
	double eucl_distance(const Vertex&, const Vertex&);
	Vertex get_turnpoint(int i){ return _tnps_tar[i];};
	void calc_leastrisk_sink();
	pair<double, vector<int>> quick_wrapup(const vector<int>& forbidden_nodes, int);
	bool check_dominance(int node, double pathrisk, double pathreward);
	void update_domilabels(int node, double pathrisk, double pathreward);
	void recursive_search(int, vector<int>, double, double, bool);
	bool is_intersected(vector<int>&, vector<int>&);
	void print_opt_sol();

};


#endif // !_PAR
