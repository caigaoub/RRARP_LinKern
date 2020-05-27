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
#include "MacroName.h"
#include <string>



using namespace std;


// --
struct BdyPoint{
	int 				_tarID; // 0,...,n, n+1
	int 				_subID; // 0,...,15
	// int 				_matID; //
};
struct EdgeY{
	BdyPoint 			_head_entry;
	BdyPoint 			_head_exit;
	BdyPoint 			_tail_entry;
	BdyPoint 			_tail_exit;
};

typedef vector<pair<BdyPoint,BdyPoint>> TOUR;


class LK{
public:
	DataHandler *												_dataset = nullptr;

	int 														_graphsize;
	vector<vector<pair<bool,double>>> * 						_G;
	vector<vector<vector<pair<double, string>>>>* 				_all_innerG;
	int 														_nb_dstzn;
	int 														_nb_targets;


	int 														_toursize;
	// TOUR														_curbesttour; // curbesttour after flip, but the entry or exit maybe not optimal 
	TOUR 														_Tour_Star; // the best known tour
	double														_Obj_Star = INF;



	void initialize(DataHandler&);
	vector<int> random_tsp_seq(int);

	void solve_linkern(TOUR& startTour, bool);
	void solve_opt2(TOUR& startTour, bool);
	
	void gain_search(TOUR& curTour, int posK, double risk_prevX, int depth, double criterion, set<int>& VisSet, double&, TOUR&, bool log_On);
	void backtrack_search(TOUR& curTour, int posK, double risk_prevX, int depth, double criterion, set<int>& VisSet, double& bestGain, TOUR& improvedTour, bool log_On);
	void greedy_search(TOUR& curTour, int posK, double risk_prevX, int depth, double criterion, set<int>& VisSet, double& bestGain, TOUR& improvedTour, bool log_On);

	void opt2_search(TOUR& curTour, int posK, double risk_prevX, int depth, double criterion, set<int>& VisSet, double& bestGain, TOUR& improvedTour, bool log_On);
	void opt3_search(TOUR& curTour, int posK, double risk_prevX, int depth, double criterion, set<int>& VisSet, double& bestGain, TOUR& improvedTour, bool log_On);

	
	vector<int> get_tsp_seq(TOUR& );
	pair<TOUR, double> solve_shortestpath_path(vector<int> & seq);
	// void solve_shortestpath_length(vector<int> & seq);

	int matID_inG(BdyPoint& bp, string brypointtype=ISENTRY);
	int tarID_inTour(int idx_tar,  TOUR & curTour);
	

	double get_Xrisk(int curtar_idx, int lattar_idx, TOUR& curtour);
	// tuple<BdyPoint, BdyPoint, double> get_Yrisk(int curtar_pos, int newngb_pos, TOUR& Tour);
	pair<EdgeY,double> EdgeY_tarK_tarJ(int posK, int posJplus1, TOUR& Tour);
	pair<EdgeY,double> EdgeY_tarKplus1_tarJplus1(int posKplus1, int posJ, TOUR& Tour);
	TOUR flip(int posK, int posJ, EdgeY& Y, EdgeY& Yc, const TOUR& curTour);
	void update_tour(TOUR&, const TOUR& tour);
	bool is_same(BdyPoint& bp1, BdyPoint& bp2);
	int locate_probe(const int& posK, const int& posJ);

	void print(TOUR& tour);
	void print_tsp(TOUR& tour);
	void print(set<int>& set1);
	void print(set<double>& set1);

	void print(BdyPoint& bp);
	void print(vector<int>& seq);

	void write_optimalpath(string filename, const TOUR& tour);

	
};













#endif 