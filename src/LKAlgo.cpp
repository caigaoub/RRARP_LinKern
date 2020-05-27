#include "LKAlgo.h"
#include "MacroName.h"
#include <algorithm>
#include <random>
#include <chrono> 
#define INF numeric_limits<double>::infinity()
void LK::initialize(DataHandler& dataset){
	this->_dataset = &dataset;
	this->_G = &(dataset._riskgraph);
	this->_all_innerG = &(dataset._all_innermatr);
	this->_nb_dstzn = dataset._nb_dstzn;
	this->_graphsize = dataset._graphsize;
	this->_toursize = 2*(dataset._nb_targets+2);
	this->_nb_targets = _dataset->_nb_targets;
} 

// int myrandom (int i) { return std::rand()%i;}
vector<int> LK::random_tsp_seq(int nb_targets){
	vector<int> seq;
	for(int i=0; i <= nb_targets+1; i++)
		seq.push_back(i);
	  std::random_device rd;
    std::mt19937 g(rd());
	std::shuffle(seq.begin()+1,seq.end()-1, g);
	// for(auto itr=tour.begin(); itr != tour.end(); itr++)
	// 	cout << *itr << " ";
	// cout << '\n';
	return seq;
}

void LK::solve_linkern(TOUR& startTour, bool log_On){
	TOUR curTour = startTour;
	double rval_X;
	int depth;
	double criterion;
	set<int> VisSet;
	double bestGain;
	TOUR improvedTour;
	int count = 0;
	string dir = ROOTDIR;
	string filename = dir+"/dat/pyFunc/optimalpath_" + to_string(count++) + ".txt";
	write_optimalpath(filename, curTour);
	// auto seqe = get_tsp_seq(startTour);
	// auto iniret = solve_shortestpath_path(seqe);
	while(true){
		depth = 0;
		criterion = 0.0;
		VisSet.clear();
		bestGain = 0.0;
		improvedTour.clear();
		if(log_On){
			cout << '\n';
			cout << "-------------------------------  Initial tour --------------------------------------" << endl;
			print_tsp(curTour);
			print(curTour);
			cout << "--------------------------------------------------------------------------------" << endl;
			cout << '\n';
		}

		int posK=0;
		for(posK = 0; posK <= _nb_targets; posK++){
			rval_X = get_Xrisk(posK, posK+1, curTour); 
			if(log_On){
				cout << '\n';
				cout << " ... FIRST BREAK X(" << tarID_inTour(posK, curTour) << "," << tarID_inTour(posK+1, curTour);
				cout << ")" << endl;
				cout << "depth: " << depth << " criterion: " << criterion << " bestGain: " << bestGain << endl;
			}
			VisSet.insert(tarID_inTour(posK, curTour));
			VisSet.insert(tarID_inTour(posK+1, curTour));
			gain_search(curTour, posK, rval_X, depth, criterion, VisSet, bestGain, improvedTour, log_On);	
			// opt2_search(curTour, posK, rval_X, depth, criterion, VisSet, bestGain, improvedTour, log_On);	

			VisSet.clear();
			if(bestGain > 0.0){
				// update curTour with improvedTour
				if(log_On){
					cout << " ---------------------------Dijkstra(b4)--------------------------------" << endl;
					print(improvedTour);
				}
				auto seq = get_tsp_seq(improvedTour);
				auto newtour = solve_shortestpath_path(seq);
				if (newtour.second < _Obj_Star){
					_Obj_Star = newtour.second;
				}
				// cout << newtour.second << ", ";
				update_tour(curTour, newtour.first); /* update curTour with newtour, which is solved by Dijkstra */
				if(true){
					// print(curTour);
					dir = ROOTDIR;
					filename = dir+"/dat/pyFunc/optimalpath_" + to_string(count++) + ".txt";
					write_optimalpath(filename, curTour);
					// cout << " ----------------------------- End -------------------------------------" << endl;
				}
				// if(count > 1)
				// 	exit(0);
				break;
			}

		}
		if(posK == _nb_targets+1)
			break;
	}
	// dijkstra();
}

/** In current tour, we consider break up X(K,K+1) and search promising gain for tarK 
    no matther whether posJ is b4 or after posK, we will add two new edges: Y(K,J) and Y(K+1, J+1).
    Y(K,J) will be added for sure if it satisfies conditions. Y(K+1,J+1) is added only as closing-up edge. If there 
    is another flip, Y(K+1,J+1) will be disconnected again. 
 */ 
void LK::gain_search(TOUR& curTour, int posK, double risk_prevX, int depth, double criterion, set<int>& VisSet, double& bestGain, TOUR& improvedTour, bool log_On){
	if(log_On){
		cout << '\n';
		cout << " ... # of flips implemented so far: " << depth << endl;
		cout << " ... "; print_tsp(curTour);
		cout << " ... visited nodes: "; print(VisSet);
		cout << " ... search gain for current probe: " << tarID_inTour(posK, curTour) << ", pos " << posK << endl;
	}

	if(depth < 2){//backtracking steps
		backtrack_search(curTour, posK, risk_prevX, depth,criterion, VisSet, bestGain, improvedTour, log_On);

	}else{
		greedy_search(curTour, posK, risk_prevX, depth, criterion, VisSet, bestGain, improvedTour, log_On);
	}
}

void LK::backtrack_search(TOUR& curTour, int posK, double risk_prevX, int depth, double criterion, set<int>& VisSet, double& bestGain, TOUR& improvedTour, bool log_On){
	if(log_On){
		cout << " ................... enter <BACKTRACK SEARCH> ..................... " << endl;
	}
	set<pair<double, int>> Jcandidates;
	double UB = criterion + risk_prevX;
	for(int posJ =0; posJ<=_nb_targets; posJ++){
		int tarJ = tarID_inTour(posJ, curTour);
		int tarJplus1 = tarID_inTour(posJ + 1, curTour);
		bool unvisited = (VisSet.find(tarJ)== VisSet.end() && VisSet.find(tarJplus1)==VisSet.end());
		if(unvisited){ // unvisited before
			auto ret_Y = EdgeY_tarK_tarJ(posK, posJ, curTour);
			if(ret_Y.second == INF){
				cerr << " distance error " << endl;
			}
			if(ret_Y.second < UB){ 
				Jcandidates.insert(make_pair(ret_Y.second, posJ));
			}		
		}
	}
	if(log_On){
		cout << " ... Admissible neighbors: ";
		for(auto itr=Jcandidates.begin(); itr != Jcandidates.end(); itr++ ){
			cout << "(" << tarID_inTour(itr->second, curTour) << "," << tarID_inTour(itr->second+1, curTour) << "), ";
		}
		cout << '\n';
	}

	if(Jcandidates.empty()){
		if(log_On){
			cout << " ... Backtrack: NO AVAILABLE NODES. RETURN" << endl;
		}
		return;
	}

	for(auto itr=Jcandidates.begin(); itr != Jcandidates.end(); itr++ ){
			int posJ = (*itr).second;
			int tarJ = tarID_inTour(posJ, curTour);
			int tarJplus1 = tarID_inTour(posJ + 1, curTour);
			bool unvisited = (VisSet.find(tarJ)== VisSet.end() && VisSet.find(tarJplus1)==VisSet.end());

			/*search promising gain by connecting K with J */
			if(log_On){
				cout << " ... <@For Loop> visited nodes: "; print(VisSet);
				cout << " *** posJ = " << posJ <<  " [break X1(" << tarID_inTour(posK,curTour) << "," << tarID_inTour(posK+1,curTour) << ")";
				cout << " and X2(" << tarID_inTour(posJ,curTour) << "," << tarID_inTour(posJ+1,curTour) << ")]";
				cout << " & [add Y(" << tarID_inTour(posK, curTour) << "," << tarID_inTour(posJ,curTour) << ")";
				cout << " and Yc(" << tarID_inTour(posK+1, curTour) << "," << tarID_inTour(posJ+1,curTour) << ")]" << endl;
			}
			auto ret_Y = EdgeY_tarK_tarJ(posK, posJ, curTour);
			if(ret_Y.second == INF){
				cerr << " distance error " << endl;
			}
			if(log_On){
				cout << " ... rval_X1 = " << risk_prevX << ", and Y(";
				cout << tarID_inTour(posK, curTour) <<"," << tarJ << ") = " << ret_Y.second << endl;
			}
			double delta = risk_prevX - ret_Y.second;
			// cout << " ... criterion is " << criterion << ", " << risk_prevX <<", " << ret_Y.second <<", " << delta << endl;
			if(criterion + delta > 0){
				// criterion += delta;
				double rval_X2 = get_Xrisk(posJ, posJ+1, curTour);
				// closing up path
				auto ret_Yc = EdgeY_tarKplus1_tarJplus1(posK+1, posJ+1, curTour); 
				if(ret_Yc.second == INF){
					cerr << " distance error " << endl;
				}
				 // make the flip 
				if(log_On){
					cout << " ... rval_X2 = " << rval_X2 << ", and Yc(";
					cout << tarID_inTour(posK+1, curTour) <<"," << tarID_inTour(posJ+1,curTour) << ") = " << ret_Yc.second << endl;
				}
				auto newTour = flip(posK, posJ, ret_Y.first, ret_Yc.first, curTour);				
				if(log_On){
					cout << " ... curTour: ";
					print(curTour);
					cout << " ... newTour: ";
					print(newTour);	
				}
				
				VisSet.insert(tarJ);
				VisSet.insert(tarJplus1);
				int posNewprobe = locate_probe(posK, posJ);
				// cout << " ===: " << criterion  << ", " << rval_X2 << ", " <<  ret_Yc.second << ", " << bestGain << endl;
				// cout << " ===: " << criterion + rval_X2 - ret_Yc.second << ", " << bestGain << endl;
				if(criterion + delta + rval_X2 - ret_Yc.second > bestGain){
					update_tour(improvedTour, newTour);
					bestGain = criterion + delta + rval_X2 - ret_Yc.second;
					if(log_On){
						cout << "*********************************************************** " << endl;
						cout << " ... !!! POSITIVE GAIN !!! Gain collected so far: " << bestGain << endl;
						cout << "*********************************************************** " << endl;
						// print(_curbesttour);	
					}
					/* search all futher flips based on current flip*/
					this->gain_search(improvedTour, posNewprobe, rval_X2, depth+1, criterion+delta, VisSet, bestGain, improvedTour, log_On);
					/* after deeper search completes, we already have improvement. So, terminate the prog and accept the first round tour 
					improvement. No need to select another pair <J,J+1> in current depth*/
					return;
				}else{
					if(log_On){
						cout << " ... current flip satisfies the criterion but not improve directly. Need search deeper " << endl;
					}
					this->gain_search(newTour, posNewprobe, rval_X2, depth+1, criterion, VisSet, bestGain, improvedTour, log_On);
					/* after deeper search is completed, we could have improvement or have nothing. If bestgain>0, we return and don't
						need to select another <J,J+1> in the same level to seek improvement; if bestgain=0, we remove current pair <J,J+1> and 
						select another pair of <J,J+1> in current depth for flipping
					*/ 
					if(bestGain > 0.0){
						return;
					}else{
						VisSet.erase(tarJ);
						VisSet.erase(tarJplus1);	
					}					
				}
			}else{
				if(log_On){
					cout << " ... negative criterion. Backtrack (depth: "<< depth << " )\n" << endl;
				}
			}
			
			/* after implementing all pairs<J,J+1> with the given FIRST BREAK,   */
			if(posJ == _nb_targets && unvisited == false){
				if(bestGain > 0.0){
					/*find a sequence of flips with positive gain*/
					if(log_On){
						cout << " ... CURRENT FIRST BREAK MAKES IMPROVEMENT. SINCE NO NODES AVIAL TO PICK, RETURN TO UPPER LEVEL!! " << endl;
						cout << '\n';
					}
					return;
				}else{
					if(log_On){
						cout << " ... NO MORE NODE AVAILABLE TO PICK AND CHOOSE ANOTHER FIRST BREAK!! " << endl;
						cout << '\n';
					}
				}
			}
		}
}

void LK::greedy_search(TOUR& curTour, int posK, double risk_prevX, int depth, double criterion, set<int>& VisSet, double& bestGain, TOUR& improvedTour, bool log_On){
	if(log_On){
		cout << " ................... enter <GREEDY SEARCH> ..................... " << endl;
	}
	set<pair<double, int>> Jcandidates;
	// int best_posJ = -1;
	double UB = criterion + risk_prevX;
	// double bestrval = INF;
	// pair<EdgeY,double> bestret;
	for(int posJ =0; posJ<=_nb_targets; posJ++){
		int tarJ = tarID_inTour(posJ, curTour);
		int tarJplus1 = tarID_inTour(posJ + 1, curTour);
		bool unvisited = (VisSet.find(tarJ)== VisSet.end() && VisSet.find(tarJplus1)==VisSet.end());
		if(unvisited){ // unvisited before
			auto ret_Y = EdgeY_tarK_tarJ(posK, posJ, curTour);
			if(ret_Y.second == INF){
				cerr << " distance error " << endl;
			}
			if(ret_Y.second < UB){ //
				// if(ret_Y.second < bestrval){
				// 	best_posJ = posJ;
				// 	bestrval = ret_Y.second;
				// 	bestret = ret_Y;
				// }
				Jcandidates.insert(make_pair(ret_Y.second, posJ));
			}		
		}
	}

	if(Jcandidates.empty()){
		if(log_On){
			cout << " ... GREEDY: NO AVAILABLE NODES. RETURN" << endl;
		}
		return;
	}else{
		// int & posJ = best_posJ;
		int posJ = (*Jcandidates.begin()).second;
		int tarJ = tarID_inTour(posJ, curTour);
		int tarJplus1 = tarID_inTour(posJ + 1, curTour);	
		/*search promising gain by connecting K with J */
		if(log_On){
			cout << " ... <@GREEDY> visited nodes: "; print(VisSet);
			cout << " *** " <<  " [break X1(" << tarID_inTour(posK,curTour) << "," << tarID_inTour(posK+1,curTour) << ")";
			cout << " and X2(" << tarID_inTour(posJ,curTour) << "," << tarID_inTour(posJ+1,curTour) << ")]";
			cout << " & [add Y(" << tarID_inTour(posK, curTour) << "," << tarID_inTour(posJ,curTour) << ")";
			cout << " and Yc(" << tarID_inTour(posK+1, curTour) << "," << tarID_inTour(posJ+1,curTour) << ")]" << endl;
		}
		// pair<EdgeY,double> & ret_Y = bestret;
		auto ret_Y = EdgeY_tarK_tarJ(posK, posJ, curTour);
		if(log_On){
			cout << " ... rval_X1 = " << risk_prevX << ", and Y(";
			cout << tarID_inTour(posK, curTour) <<"," << tarJ << ") = " << ret_Y.second << endl;
		}
		double delta = risk_prevX - ret_Y.second;
		// cout << " ... criterion is " << criterion + delta << endl;
		criterion += delta;
		double rval_X2 = get_Xrisk(posJ, posJ+1, curTour);
		// closing up path
		auto ret_Yc = EdgeY_tarKplus1_tarJplus1(posK+1, posJ+1, curTour); 
		if(ret_Yc.second == INF){
			cerr << " distance error " << endl;
		}
		 // make the flip 
		if(log_On){
			cout << " ... rval_X2 = " << rval_X2 << ", and Yc(";
			cout << tarID_inTour(posK+1, curTour) <<"," << tarID_inTour(posJ+1,curTour) << ") = " << ret_Yc.second << endl;
		}
		auto newTour = flip(posK, posJ, ret_Y.first, ret_Yc.first, curTour);				
		if(log_On){
			cout << " ... curTour: ";
			print(curTour);
			cout << " ... newTour: ";
			print(newTour);	
		}
		
		VisSet.insert(tarJ);
		VisSet.insert(tarJplus1);
		int posNewprobe = locate_probe(posK, posJ);
		if(criterion + rval_X2 - ret_Yc.second > bestGain){
			update_tour(improvedTour, newTour);
			bestGain = criterion + rval_X2 - ret_Yc.second;
			if(log_On){
				cout << "*********************************************************** " << endl;
				cout << " ... !!! POSITIVE GAIN !!! Gain collected so far: " << bestGain << endl;
				cout << "*********************************************************** " << endl;
			}
			this->gain_search(improvedTour, posNewprobe, rval_X2, depth+1, criterion, VisSet, bestGain, improvedTour, log_On);
			if(log_On){
				cout << " ... find improvement " << endl;
			}	
			// return;
		}else{
			if(log_On){
				cout << " ... current flip satisfies the criterion but not improve directly. Need search deeper " << endl;
			}
			this->gain_search(newTour, posNewprobe, rval_X2, depth+1, criterion, VisSet, bestGain, improvedTour, log_On);
			if(bestGain > 0.0){
				if(log_On){
					cout << " ... find improvement " << endl;
				}	
				// return;
			}else{
				// VisSet.erase(tarJ);
				// VisSet.erase(tarJplus1);
				if(log_On){
					cout << " ... no improvement still" << endl;
				}	
			}					
		}
	}

}



void LK::solve_opt2(TOUR& startTour, bool log_On){
	TOUR curTour = startTour;
	double rval_X;
	int depth;
	double criterion;
	set<int> VisSet;
	double bestGain;
	TOUR improvedTour;
	int count = 0;
	string dir = ROOTDIR;
	string filename = dir+"/dat/pyFunc/optimalpath_" + to_string(count++) + ".txt";
	// write_optimalpath(filename, curTour);
	// auto seqe = get_tsp_seq(startTour);
	// auto iniret = solve_shortestpath_path(seqe);
	int nb_improvements = 0;
	while(true){

		if(log_On){
			cout << '\n';
			cout << "-------------------------------  Initial tour --------------------------------------" << endl;
			print_tsp(curTour);
			print(curTour);
			cout << "--------------------------------------------------------------------------------" << endl;
			cout << '\n';
		}

		int posK=0;
		for(posK = 0; posK <= _nb_targets; posK++){
			depth = 0;
			criterion = 0.0;
			VisSet.clear();
			bestGain = 0.0;
			improvedTour.clear();
			rval_X = get_Xrisk(posK, posK+1, curTour); 
			if(log_On){
				cout << '\n';
				cout << " ... FIRST BREAK X(" << tarID_inTour(posK, curTour) << "," << tarID_inTour(posK+1, curTour);
				cout << ")" << endl;
				cout << "depth: " << depth << " criterion: " << criterion << " bestGain: " << bestGain << endl;
			}
			VisSet.insert(tarID_inTour(posK, curTour));
			VisSet.insert(tarID_inTour(posK+1, curTour));
			opt2_search(curTour, posK, rval_X, depth, criterion, VisSet, bestGain, improvedTour, log_On);	

			VisSet.clear();
			if(bestGain > 0.0){
				nb_improvements ++;
				if(nb_improvements > 10){
					nb_improvements = 0;
					auto seq = get_tsp_seq(improvedTour);
					auto newtour = solve_shortestpath_path(seq);
					update_tour(curTour, newtour.first); /* update curTour with newtour, which is solved by Dijkstra */
					if (newtour.second < _Obj_Star){
						_Obj_Star = newtour.second;
					}
				}else{
					update_tour(curTour, improvedTour); 
				}

				if(true){
					// print(curTour);
					dir = ROOTDIR;
					filename = dir+"/dat/pyFunc/optimalpath_" + to_string(count++) + ".txt";
					// write_optimalpath(filename, curTour);
					// cout << " ----------------------------- End -------------------------------------" << endl;
				}
				// if(count > 1)
				// 	exit(0);
				break;
			}

		}
		if(improvedTour.size() != 0){
			auto seq = get_tsp_seq(improvedTour);
			auto newtour = solve_shortestpath_path(seq);
			_Obj_Star = newtour.second;
		}
		if(posK == _nb_targets+1)
			break;
	}
	// dijkstra();
}
void LK::opt2_search(TOUR& curTour, int posK, double risk_prevX, int depth, double criterion, set<int>& VisSet, double& bestGain, TOUR& improvedTour, bool log_On){
	if(log_On){
		cout << '\n';
		cout << " ... # of flips implemented so far: " << depth << endl;
		cout << " ... "; print_tsp(curTour);
		cout << " ... visited nodes: "; print(VisSet);
		cout << " ... search gain for current probe: " << tarID_inTour(posK, curTour) << ", pos " << posK << endl;
	}

	set<pair<double, int>> Jcandidates;
	double UB = criterion + risk_prevX;
	for(int posJ =0; posJ<=_nb_targets; posJ++){
		int tarJ = tarID_inTour(posJ, curTour);
		int tarJplus1 = tarID_inTour(posJ + 1, curTour);
		bool unvisited = (VisSet.find(tarJ)== VisSet.end() && VisSet.find(tarJplus1)==VisSet.end());
		if(unvisited){ // unvisited before
			auto ret_Y = EdgeY_tarK_tarJ(posK, posJ, curTour);
			if(ret_Y.second == INF){
				cerr << " distance error " << endl;
			}
			if(ret_Y.second < UB){ 
				Jcandidates.insert(make_pair(ret_Y.second, posJ));
			}		
		}
	}

	if(Jcandidates.empty()){
		if(log_On){
			cout << " ... Backtrack: NO AVAILABLE NODES. RETURN" << endl;
		}
		return;
	}

	for(auto itr=Jcandidates.begin(); itr != Jcandidates.end(); itr++ ){
			int posJ = (*itr).second;
			int tarJ = tarID_inTour(posJ, curTour);
			int tarJplus1 = tarID_inTour(posJ + 1, curTour);
			bool unvisited = (VisSet.find(tarJ)== VisSet.end() && VisSet.find(tarJplus1)==VisSet.end());

			/*search promising gain by connecting K with J */
			if(log_On){
				cout << " ... <@For Loop> visited nodes: "; print(VisSet);
				cout << " *** posJ = " << posJ <<  " [break X1(" << tarID_inTour(posK,curTour) << "," << tarID_inTour(posK+1,curTour) << ")";
				cout << " and X2(" << tarID_inTour(posJ,curTour) << "," << tarID_inTour(posJ+1,curTour) << ")]";
				cout << " & [add Y(" << tarID_inTour(posK, curTour) << "," << tarID_inTour(posJ,curTour) << ")";
				cout << " and Yc(" << tarID_inTour(posK+1, curTour) << "," << tarID_inTour(posJ+1,curTour) << ")]" << endl;
			}
			auto ret_Y = EdgeY_tarK_tarJ(posK, posJ, curTour);
			if(ret_Y.second == INF){
				cerr << " distance error " << endl;
			}
			if(log_On){
				cout << " ... rval_X1 = " << risk_prevX << ", and Y(";
				cout << tarID_inTour(posK, curTour) <<"," << tarJ << ") = " << ret_Y.second << endl;
			}
			double delta = risk_prevX - ret_Y.second;
			// if(delta > 0){
				// criterion += delta;
				double rval_X2 = get_Xrisk(posJ, posJ+1, curTour);
				// closing up path
				auto ret_Yc = EdgeY_tarKplus1_tarJplus1(posK+1, posJ+1, curTour); 
				if(ret_Yc.second == INF){
					cerr << " distance error " << endl;
				}
				 // make the flip 
				if(log_On){
					cout << " ... rval_X2 = " << rval_X2 << ", and Yc(";
					cout << tarID_inTour(posK+1, curTour) <<"," << tarID_inTour(posJ+1,curTour) << ") = " << ret_Yc.second << endl;
				}
				auto newTour = flip(posK, posJ, ret_Y.first, ret_Yc.first, curTour);				
				if(log_On){
					cout << " ... curTour: ";
					print(curTour);
					cout << " ... newTour: ";
					print(newTour);	
				}
				
				// VisSet.insert(tarJ);
				// VisSet.insert(tarJplus1);
				// int posNewprobe = locate_probe(posK, posJ);
				// cout << " ===: " << criterion  << ", " << rval_X2 << ", " <<  ret_Yc.second << ", " << bestGain << endl;
				// cout << " ===: " << criterion + rval_X2 - ret_Yc.second << ", " << bestGain << endl;
				if( delta + rval_X2 - ret_Yc.second > bestGain){
					update_tour(improvedTour, newTour);
					bestGain = criterion + delta + rval_X2 - ret_Yc.second;
					if(log_On){
						cout << "*********************************************************** " << endl;
						cout << " ... !!! POSITIVE GAIN !!! Gain collected so far: " << bestGain << endl;
						cout << "*********************************************************** " << endl;
						// print(_curbesttour);	
					}
					/* search all futher flips based on current flip*/
					// this->gain_search(improvedTour, posNewprobe, rval_X2, depth+1, criterion+delta, VisSet, bestGain, improvedTour, log_On);
					/* after deeper search completes, we already have improvement. So, terminate the prog and accept the first round tour 
					improvement. No need to select another pair <J,J+1> in current depth*/
					return;
				}
			// }
			
			
		}
}



void LK::opt3_search(TOUR& curTour, int posK, double risk_prevX, int depth, double criterion, set<int>& VisSet, double& bestGain, TOUR& improvedTour, bool log_On){
	if(log_On){
		cout << '\n';
		cout << " ... # of flips implemented so far: " << depth << endl;
		cout << " ... "; print_tsp(curTour);
		cout << " ... visited nodes: "; print(VisSet);
		cout << " ... search gain for current probe: " << tarID_inTour(posK, curTour) << ", pos " << posK << endl;
	}

	if(depth < 3){//backtracking steps
		backtrack_search(curTour, posK, risk_prevX, depth,criterion, VisSet, bestGain, improvedTour, log_On);
	}
}

int LK::matID_inG(BdyPoint& bp, string brypointtype){
	if(bp._tarID == 0){
		return 0;
	}else if(bp._tarID == _nb_targets+1){
		return _graphsize-1;
	}else{
		if(brypointtype.compare(ISEXIT)==0){
			return (bp._tarID-1)*2*_nb_dstzn + _nb_dstzn + 1 + bp._subID;		
		}
		return (bp._tarID-1)*2*_nb_dstzn + 1 + bp._subID;	
	}
}


int LK::tarID_inTour(int idx_tar,  TOUR & curTour){
	return curTour[idx_tar].first._tarID;
}


TOUR LK::flip(int posK, int posJ, EdgeY& Y, EdgeY& Yc, const TOUR& curTour){
	TOUR newTour;
	newTour.resize(_nb_targets+2);
	if(posK < posJ){ // current tour: 0 ~~~>K--K+1~~~~>J--J+1~~~~~>n+1
		for(int i=0; i<=posK; i++){
			newTour[i] = curTour[i];
		}
		for(int i=posK+1; i<=posJ; i++){ // revert [K+1,...,J]
			newTour[i] = make_pair(curTour[posJ+posK+1-i].second, curTour[posJ+posK+1-i].first);
		}
		for(int i=posJ+1; i<= _nb_targets+1; i++){
			newTour[i] = curTour[i];
		}
		/* follow the sequence of Y: head_entry@K -> head_exit@K -> tail_entry@J -> tail_exit@J */
		newTour[posK] = {Y._head_entry, Y._head_exit};
		newTour[posK+1] = {Y._tail_entry, Y._tail_exit};
		/* follow the sequence of Yc: head_entry@K -> head_exit@K -> tail_entry@J -> tail_exit@J */
		newTour[posJ] = {Yc._head_entry, Yc._head_exit};
		newTour[posJ+1] = {Yc._tail_entry, Yc._tail_exit};

	}
	if(posJ < posK){// current tour: 0 ~~~>J--J+1~~~~>K--K+1~~~~~>n+1
		for(int i=0; i<=posJ; i++){
			newTour[i] = curTour[i];
		}
		for(int i=posJ+1; i<=posK; i++){// revert [J+1,...,K]
			newTour[i] = make_pair(curTour[posK+posJ+1-i].second, curTour[posK+posJ+1-i].first);
		}
		for(int i=posK+1; i<= _nb_targets+1; i++){
			newTour[i] = curTour[i];
		}
		/*reverse of Y: head_entry@K <- head_exit@K <- tail_entry@J <- tail_exit@J  */
		newTour[posJ] = {Y._tail_exit, Y._tail_entry};
		newTour[posJ+1] = {Y._head_exit, Y._head_entry}; //careful
		/*reverse of Yc: head_entry@K <- head_exit@K <- tail_entry@J <- tail_exit@J  */
		newTour[posK] = {Yc._tail_exit, Yc._tail_entry};
		newTour[posK+1] = { Yc._head_exit, Yc._head_entry};
	}
	return newTour;
}

/* EdgeY will be like: ak -> bk ~~> aj ->bj (if posK < posJ) 
		ak: current entry on tarK;
		bk: the new exit on tarK;
		aj: the new entry on tarJ;
		bj: current entry on tarJ. (now it acts as tarJ's new exit)
	  or aj -> bj ~~> ak ->bk (if posJ < posK)
	  	aj: current entry on tarJ;
		bj: the new exit on tarJ;
		ak: the new entry on tarK;
		bk: current entry on tarK. (now it acts as tarK's new exit)
	Thus, pick two enties as endpoints of EdgeY
	*/
pair<EdgeY,double> LK::EdgeY_tarK_tarJ(int posK, int posJ, TOUR& Tour){

	// cout << " ... FUNCTION: add Y" << endl;
	int i_star = -1, j_star = -1;
	BdyPoint bpi, bpj;
	double riskval = INF, rval_tmp = 0;
	/* 	the two endpoints of new EdgeY(K,J) are fixed at first. Head endpoint is the entry of tarK 
		and the tail endpoint is also the entry of the tarJ (due to the reverse) */
	BdyPoint& head_Y = Tour[posK].first;
	BdyPoint& tail_Y = Tour[posJ].first;
	// cout << " ...... HEAD "; print(head_Y);
	// cout << " ...... TAIL "; print(tail_Y);
	if(posK == 0){ // if K is start depot
		for(int j=0; j<_nb_dstzn; j++){//on tarJ
			bpj = {tarID_inTour(posJ, Tour), j};
			if(is_same(bpj, tail_Y))
				continue;
			rval_tmp = (*_G)[0][matID_inG(bpj)].second; // outer trajc
			rval_tmp += (*_G)[matID_inG(bpj)][matID_inG(tail_Y, ISEXIT)].second; // inner trajc
			// cout << j << " r1=" << (*_G)[0][matID_inG(bpj)].second << " r2=" << (*_G)[matID_inG(bpj)][matID_inG(tail_Y,ISEXIT)].second << endl;
			if(rval_tmp < riskval){
				riskval = rval_tmp;
				j_star = j;
			}
		}
		// cout << " ...... j*_on_tarJ=" << j_star << " rval=" << riskval << endl;
		bpj = {Tour[posJ].first._tarID, j_star};
		EdgeY Y = {head_Y, head_Y, bpj, tail_Y};
		return make_pair(Y, riskval);		
	}else if(posJ == 0){	
		for(int i=0; i<_nb_dstzn; i++){ //on tarK
			bpi = {tarID_inTour(posK, Tour), i};
			if(is_same(bpi, head_Y))
				continue;
			rval_tmp = (*_G)[matID_inG(head_Y)][matID_inG(bpi, ISEXIT)].second; // outer
			rval_tmp += (*_G)[matID_inG(bpi)][0].second; // inner 
			// cout << i << " r1=" << (*_G)[matID_inG(head_Y)][matID_inG(bpi, ISEXIT)].second << " r2=" << (*_G)[matID_inG(bpi)][0].second << endl;
			
			if(rval_tmp < riskval){
				riskval = rval_tmp;
				i_star = i;
			}
		}
		// cout << " ...... i*_on_tarK =" << i_star << " rval=" << riskval << endl;
		bpi = {Tour[posK].first._tarID, i_star};
		// bpj = {Tour[posJ].first._tarID, j_star};
		EdgeY Y = {head_Y, bpi, tail_Y, tail_Y};
		return make_pair(Y, riskval);
	}else{
		for(int i=0; i<_nb_dstzn; i++){ //on tarK
			bpi = {tarID_inTour(posK, Tour), i};
			if(is_same(bpi, head_Y))
				continue;
			for(int j=0; j<_nb_dstzn; j++){//on tarJplus1
				bpj = {tarID_inTour(posJ, Tour), j};
				if(is_same(bpj, tail_Y))
					continue;
				rval_tmp = (*_G)[matID_inG(head_Y)][matID_inG(bpi, ISEXIT)].second;
				rval_tmp += (*_G)[matID_inG(bpi, ISEXIT)][matID_inG(bpj)].second;
				rval_tmp += (*_G)[matID_inG(bpj)][matID_inG(tail_Y,ISEXIT)].second;
				// cout << i << "-"<< j << " r1=" << (*_G)[matID_inG(head_Y)][matID_inG(bpi, ISEXIT)].second << " r2=" << (*_G)[matID_inG(bpi, ISEXIT)][matID_inG(bpj)].second;
				// cout << " r3=" << (*_G)[matID_inG(bpj)][matID_inG(tail_Y,ISEXIT)].second << endl;
				if(rval_tmp < riskval){
					riskval = rval_tmp;
					i_star = i;
					j_star = j;
				}
			}
		}
		// cout << " ...... i*_on_tarK=" << i_star << " j*_on_tarJ=" << j_star << " rval=" << riskval << endl;
		bpi = {Tour[posK].first._tarID, i_star};
		bpj = {Tour[posJ].first._tarID, j_star};
		EdgeY Y = {head_Y, bpi, bpj, tail_Y};
		return make_pair(Y, riskval);
	}

	cerr << "should not reach here !!!!!!!!!!" << endl;
}

pair<EdgeY,double> LK::EdgeY_tarKplus1_tarJplus1(int posKplus1, int posJplus1, TOUR& Tour){
	int i_star = -1, j_star = -1;
	BdyPoint bpi, bpj;
	double riskval = INF, rval_tmp = 0;
	BdyPoint& head_Yc = Tour[posKplus1].second;
	BdyPoint& tail_Yc = Tour[posJplus1].second;
	// cout << " ...... HEAD Yc "; print(head_Yc);
	// cout << " ...... TAIL Yc "; print(tail_Yc);
	if(posKplus1 == _nb_targets+1){ // if K+1 is end depot
		for(int j=0; j<_nb_dstzn; j++){//on tarJ
			bpj = {tarID_inTour(posJplus1, Tour), j};
			if(is_same(bpj, tail_Yc))
				continue;
			rval_tmp = (*_G)[_graphsize-1][matID_inG(bpj, ISEXIT)].second; // outer
			rval_tmp += (*_G)[matID_inG(bpj, ISEXIT)][matID_inG(tail_Yc)].second;
			// cout << j << " ...... r1=" << (*_G)[_graphsize-1][matID_inG(bpj,ISEXIT)].second << " r2=" << (*_G)[matID_inG(bpj,ISEXIT)][matID_inG(tail_Yc)].second << endl;
			if(rval_tmp < riskval){
				riskval = rval_tmp;
				j_star = j;
			}
		}
		// cout << " ...... j*_on_tarJ=" << j_star << ", rval=" << riskval << endl;
		bpj = {Tour[posJplus1].first._tarID, j_star};
		EdgeY Y = {head_Yc, head_Yc, bpj, tail_Yc};
		return make_pair(Y, riskval);		
	}else if(posJplus1 == _nb_targets+1){	
		for(int i=0; i<_nb_dstzn; i++){ //on tarK
			bpi = {tarID_inTour(posKplus1, Tour), i};
			if(is_same(bpi, head_Yc))
				continue;
			rval_tmp = (*_G)[matID_inG(head_Yc)][matID_inG(bpi, ISEXIT)].second;
			rval_tmp += (*_G)[matID_inG(bpi, ISEXIT)][_graphsize-1].second;
			// cout << i << " r1=" << (*_G)[matID_inG(head_Y)][matID_inG(bpi, ISEXIT)].second << " r2=" << (*_G)[matID_inG(bpi)][0].second << endl;
			if(rval_tmp < riskval){
				riskval = rval_tmp;
				i_star = i;
			}
		}
		// cout << " ...... i*_on_tarK =" << i_star << " rval=" << riskval << endl;
		bpi = {Tour[posKplus1].first._tarID, i_star};
		// bpj = {Tour[posJ].first._tarID, j_star};
		EdgeY Yc = {head_Yc, bpi, tail_Yc, tail_Yc};
		return make_pair(Yc, riskval);
	}else{
		for(int i=0; i<_nb_dstzn; i++){ //on tarK
			bpi = {tarID_inTour(posKplus1, Tour), i};
			if(is_same(bpi, head_Yc))
				continue;
			for(int j=0; j<_nb_dstzn; j++){//on tarJplus1
				bpj = {tarID_inTour(posJplus1, Tour), j};
				if(is_same(bpj, tail_Yc))
					continue;
				rval_tmp = (*_G)[matID_inG(head_Yc)][matID_inG(bpi, ISEXIT)].second;
				rval_tmp += (*_G)[matID_inG(bpi, ISEXIT)][matID_inG(bpj)].second;
				rval_tmp += (*_G)[matID_inG(bpj)][matID_inG(tail_Yc,ISEXIT)].second;
				// cout << i << "-"<< j << " r1=" << (*_G)[matID_inG(head_Y)][matID_inG(bpi, ISEXIT)].second << " r2=" << (*_G)[matID_inG(bpi, ISEXIT)][matID_inG(bpj)].second;
				// cout << " r3=" << (*_G)[matID_inG(bpj)][matID_inG(tail_Y,ISEXIT)].second << endl;
				if(rval_tmp < riskval){
					riskval = rval_tmp;
					i_star = i;
					j_star = j;
				}
			}
		}
		// cout << " ...... i*_on_tarK=" << i_star << " j*_on_tarJ=" << j_star << " rval=" << riskval << endl;
		bpi = {Tour[posKplus1].first._tarID, i_star};
		bpj = {Tour[posJplus1].first._tarID, j_star};
		EdgeY Yc = {head_Yc, bpi, bpj, tail_Yc};
		return make_pair(Yc, riskval);
	}

	cerr << "should not reach here !!!!!!!!!!" << endl;
}

double LK::get_Xrisk(int curtar_idx, int lattar_idx, TOUR& curtour){
	// int curtar_ID = curtour[curtar_idx].first._tarID;
	// int lattar_ID = curtour[lattar_idx].first._tarID;
	double riskval = INF;
	if(curtar_idx == 0){
		int entry_lat = matID_inG(curtour[lattar_idx].first);
		int exit_lat = matID_inG(curtour[lattar_idx].second, ISEXIT);
		riskval  = (*_G)[0][entry_lat].second + (*_G)[entry_lat][exit_lat].second;
	}else if(lattar_idx == _nb_targets+1 ){
		int entry_cur = matID_inG(curtour[curtar_idx].first);
		int exit_cur = matID_inG(curtour[curtar_idx].second, ISEXIT);
		riskval = (*_G)[entry_cur][exit_cur].second + (*_G)[exit_cur][_graphsize-1].second;
	}else{
		int entry_cur = matID_inG(curtour[curtar_idx].first);
		int exit_cur = matID_inG(curtour[curtar_idx].second, ISEXIT);
		int entry_lat = matID_inG(curtour[lattar_idx].first);
		int exit_lat = matID_inG(curtour[lattar_idx].second, ISEXIT);
		riskval = (*_G)[entry_cur][exit_cur].second + (*_G)[exit_cur][entry_lat].second + (*_G)[entry_lat][exit_lat].second;
	}
	// cout << " ... X(" << curtour[curtar_idx].first._tarID << "," << curtour[lattar_idx].first._tarID << "): "<< riskval << endl;
	return riskval;
}

int LK::locate_probe(const int& posK, const int& posJ){
	if(posK ==posJ)
		cerr << "wrong position pick, (line286)! " << endl; 
	if(posK < posJ){
		return posJ;
	}else{
		return posK;
	}
}

bool LK::is_same(BdyPoint& bp1, BdyPoint& bp2){
	if(bp1._tarID == bp2._tarID && bp1._subID == bp2._subID)
		return true;
	else
		return false;
}

void LK::update_tour(TOUR& improvedTour, const TOUR& tour){
	if(improvedTour.empty())
		improvedTour.resize(_nb_targets+2);
	for(int i=0; i<= _nb_targets+1; i++){
		improvedTour[i] = tour[i];
	}
}


void LK::print(TOUR& tour){
	// cout << "Tour: ";
	for(auto itr=tour.begin(); itr != tour.end(); itr++){
		if(itr == tour.end()-1){
			cout << "["<<itr->first._tarID << ',' << itr->first._subID << "->" << itr->second._tarID << ',' << itr->second._subID << "]";
		}else{
			cout << "["<<itr->first._tarID << ',' << itr->first._subID << "->" << itr->second._tarID << ',' << itr->second._subID << "]" << "==>";
		}
	}
	cout << '\n';
}

void LK::print_tsp(TOUR& tour){
	// cout << "tsp: ";
	for(auto itr=tour.begin(); itr != tour.end(); itr++){
		cout << itr->first._tarID << ' ';
	}
	cout << '\n';
}

void LK::print(vector<int>& seq){
	// cout << "tsp: ";
	for(auto itr=seq.begin(); itr != seq.end(); itr++){
		cout << *itr << ' ';
	}
	cout << '\n';
}

void LK::print(set<int>& set1){
	if(set1.empty()){
		cout<< "the set is empty !!" << endl;
	}else{
		for(auto itr=set1.begin(); itr != set1.end(); itr++)
			cout << *itr << ' ';
		cout << endl;
	}

}

void LK::print(set<double>& set1){
	if(set1.empty()){
		cout<< "the set is empty !!" << endl;
	}else{
		for(auto itr=set1.begin(); itr != set1.end(); itr++)
			cout << *itr << ' ';
		cout << endl;
	}

}

void LK::print(BdyPoint& bp){
	cout << "BdyPoint: " << bp._tarID << "," << bp._subID << " ";
	_dataset->_points[bp._tarID][bp._subID].print();
}

vector<int> LK::get_tsp_seq(TOUR& tour){
	vector<int> seq;
	for(int i=0; i<=_nb_targets+1; i++){
		seq.push_back(tarID_inTour(i, tour));
	}
	return seq;
}


pair<TOUR, double> LK::solve_shortestpath_path(vector<int> & seq) {
	int idxmat_1, idxmat_2, idxmat_3;
	vector<double> entry_dist(_nb_dstzn, INF);
	vector<double> exit_dist(_nb_dstzn, INF);
	vector<double> INF_dist(_nb_dstzn, INF);
	double dist_endDepot = INF;
	double dist = 0.0;
	vector<vector<BdyPoint>> prevnodes_tars(2 * _dataset->_nb_targets); //pair<target, turning point>
	for(int i = 0;i< 2* _nb_targets; i++){
		prevnodes_tars[i].resize(_nb_dstzn);
	}

	vector<vector<double>> SDS;
	SDS.resize(_dataset->_nb_targets+2);
	BdyPoint prevnode_arrival;
	int col_layer = 0;
	// i) SDS[start depot] = 0.0
	SDS[0].resize(1);
	SDS[0][0] = 0.0;
	// ii) SDS[*] of boundary points of the first visited target
	int pos = 1;
	int idx_circle = seq[pos];
	SDS[1].resize(_nb_dstzn);
	idxmat_1 = (idx_circle - 1) * 2 * _nb_dstzn + 1; //first entry of the first visited target
	for (int i = 0; i < _nb_dstzn; i++)	{
		// cout << i << " = " << (*_G)[0][idxmat_1 + i].second << endl;
		if ((*_G)[0][idxmat_1 + i].first == true) {
			entry_dist[i] = (*_G)[0][idxmat_1 + i].second;
			SDS[1][i] = entry_dist[i];
			prevnodes_tars[col_layer][i] = {col_layer,0};
		}
	}

	col_layer++;
	idxmat_2 = (idx_circle - 1) * 2 * _nb_dstzn + _nb_dstzn + 1;
	for (int i = 0; i < _nb_dstzn; i++) {
		for (int j = 0; j < _nb_dstzn; j++) {
			if ((*_G)[idxmat_1 + i][idxmat_2 + j].first == true) {
				dist = entry_dist[i] + (*_G)[idxmat_1 + i][idxmat_2 + j].second;
				if (dist < exit_dist[j]){
					exit_dist[j] = dist;
					prevnodes_tars[col_layer][j] = {idx_circle,i};
				}
			}
		}
	}
	entry_dist = INF_dist;
	// iii) nodes from seq[2] to seq[num_circles]
	int idx_fr_circle = 0;
	int idx_lat_circle = 0;
	col_layer++;

	for(pos = 2; pos <= _dataset->_nb_targets; pos++) {
		SDS[pos].resize(_nb_dstzn);
		idx_fr_circle = seq[pos - 1]; 
		idx_lat_circle = seq[pos]; 
		idxmat_1 = (idx_fr_circle - 1) * 2 * _nb_dstzn + _nb_dstzn + 1;
		idxmat_2 = (idx_lat_circle - 1) * 2 * _nb_dstzn + 1;

		for (int i = 0; i < _nb_dstzn; i++) {
			for (int j = 0; j < _nb_dstzn; j++) {
				if ((*_G)[idxmat_1 + i][idxmat_2 + j].first == true) {
					dist = exit_dist[i] + (*_G)[idxmat_1 + i][idxmat_2 + j].second;
					if (dist < entry_dist[j]){
						entry_dist[j] = dist;
						prevnodes_tars[col_layer][j] = {idx_fr_circle, i};
					}
				}
			}
		}
		col_layer++;

		for (int i = 0; i < _nb_dstzn; i++) {
			SDS[pos][i] = entry_dist[i];
		}
		exit_dist = INF_dist;
		idxmat_3 = (idx_lat_circle - 1) * 2 * _nb_dstzn + _nb_dstzn + 1;
		// update exits
		for (int i = 0; i < _nb_dstzn; i++) {
			for (int j = 0; j < _nb_dstzn; j++) {
				if ((*_G)[idxmat_2 + i][idxmat_3 + j].first == true) {
					dist = entry_dist[i] + (*_G)[idxmat_2 + i][idxmat_3 + j].second;
					if (dist < exit_dist[j]){
						exit_dist[j] = dist;
						prevnodes_tars[col_layer][j] = {idx_lat_circle, i};
					}
				}
			}
		}
		col_layer++;
		entry_dist = INF_dist;

	}

	// iv) the last target <-> the end depot
	idx_circle = seq[_dataset->_nb_targets];
	SDS[_dataset->_nb_targets + 1].resize(1);
	idxmat_1 = (idx_circle - 1) * 2 * _nb_dstzn + _nb_dstzn + 1;
	// cout << idx_circle << ", " << _graphsize << endl;
	for (int i = 0; i < _nb_dstzn; i++) {
		if ((*_G)[idxmat_1 + i][_graphsize - 1].first == true) {
			dist = exit_dist[i] + (*_G)[idxmat_1 + i][_graphsize - 1].second;
			if (dist < dist_endDepot){
				dist_endDepot = dist;
				prevnode_arrival = {idx_circle, i};
			}
		}
	}

	SDS[_dataset->_nb_targets + 1][0] = dist_endDepot;

	deque<pair<BdyPoint,Vertex>> OptPath;
	BdyPoint bp = {_nb_targets+1, 0};
	OptPath.push_front(make_pair(bp, _dataset->_points[_dataset->_nb_targets+1][0]));
	OptPath.push_front(make_pair(prevnode_arrival, _dataset->_points[prevnode_arrival._tarID][prevnode_arrival._subID]));
	// int pre_tar = prevnode_arrival.first;
	int tpidx_pre_tar = prevnode_arrival._subID;
	for(int i = 2*_nb_targets-1; i >= 0; i--){
		OptPath.push_front(make_pair(prevnodes_tars[i][tpidx_pre_tar], _dataset->_points[prevnodes_tars[i][tpidx_pre_tar]._tarID][prevnodes_tars[i][tpidx_pre_tar]._subID]));
		tpidx_pre_tar = prevnodes_tars[i][tpidx_pre_tar]._subID;
	}

	// deque --> vector
	TOUR tour_star(_dataset->_nb_targets+2);
	bp = {0, 0};
	tour_star[0] = make_pair(bp, bp);
	auto itr = OptPath.begin()+1;
	for(int t=1; t<= _nb_targets; t++){
		tour_star[t] = make_pair(itr->first,(itr+1)->first);
		itr += 2;
	}
	bp = {_nb_targets+1,0};
	tour_star[_nb_targets+1] = make_pair(bp, bp);

	
	if(false){
		cout << " =====>> optimal entries and exits selected: " << '\n';
		for(auto itr = OptPath.begin(); itr != OptPath.end(); itr++){
			cout << " ... " << (*itr).first._tarID << ": (" << (*itr).second._x << ',' << (*itr).second._y << ")\n";
		}
		// string dir = ROOTDIR;
		// ofstream file_OptPath(dir+"/dat/pyFunc/optimalpath.txt");
		// for(auto itr = OptPath.begin(); itr != OptPath.end(); itr++){
		// 		file_OptPath << (*itr).first._tarID << '\t' << (*itr).second._x << '\t' << (*itr).second._y << '\n';
		// }
		// for(auto itr = OptPath.begin()+1; itr != OptPath.end()-1; itr+=2){
		// 	int tar = (*itr).first._tarID;
		// 	// cout << tar << " = " << (*itr).first._subID << ", " << (*(itr+1)).first._subID << " --> ";
		// 	// cout << (*_all_innerG)[tar-1][(*itr).first._subID][(*(itr+1)).first._subID].second << endl;
		// 	file_OptPath << (*_all_innerG)[tar-1][(*itr).first._subID][(*(itr+1)).first._subID].second << '\n';
		// }
		// file_OptPath.close();
	}
	return make_pair(tour_star, SDS[_nb_targets+1][0]);
}

void LK::write_optimalpath(string filename, const TOUR& tour){
	ofstream file_OptPath(filename);
	auto itr = tour.begin();
	int tar  = (*itr).first._tarID;
	file_OptPath << tar << '\t' << _dataset->_points[tar][(*itr).first._subID]._x << '\t' << _dataset->_points[tar][(*itr).first._subID]._y << '\n';
	itr++;
	for(; itr != tour.end()-1;){
		tar  = (*itr).first._tarID;
		file_OptPath << tar << '\t' << _dataset->_points[tar][(*itr).first._subID]._x << '\t' << _dataset->_points[tar][(*itr).first._subID]._y << '\n';
		file_OptPath << tar << '\t' << _dataset->_points[tar][(*itr).second._subID]._x << '\t' << _dataset->_points[tar][(*itr).second._subID]._y << '\n';
		itr++;
	}
	tar  = (*itr).first._tarID;
	file_OptPath << _nb_targets+1 << '\t' << _dataset->_points[tar][(*itr).first._subID]._x << '\t' << _dataset->_points[tar][(*itr).first._subID]._y << '\n';
	tar = 0;
	for(itr = tour.begin()+1; itr != tour.end()-1; itr++){
		tar = (*itr).first._tarID;
		int e1 = (*itr).first._subID;
		int e2 = (*itr).second._subID;
		// file_OptPath << tar << '\t' << (*_all_innerG)[tar-1][e1][e2].second << '\n';
		file_OptPath << (*_all_innerG)[tar-1][e1][e2].second << '\n';

		tar++;
	}
	file_OptPath.close();
}

