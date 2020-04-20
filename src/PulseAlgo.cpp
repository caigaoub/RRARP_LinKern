#include "PulseAlgo.h"
#define _USE_MATH_DEFINES
#include<bits/stdc++.h>

void Pulse::set_parameters(const int source, const int sink, const double demand_pct){
	this->_source = source;
	this->_sink = sink;
	this->_demand_pct = demand_pct;
	// this->_demand = this->_total_reward * this->_demand_pct;
}

void Pulse::initialize_triGraph(DataHandler& dataset, int tar, int lmt_obps, string graphfile){
	this->_dataset = &dataset;
	this->_taridx = tar;
	this->_graphsize = lmt_obps;

	_rewards_etd.resize(_graphsize, 0);
	_total_reward = 0.0;
	for(int i=0; i< _graphsize; i++){
		_rewards_etd[i] = _dataset->_all_obps[_taridx][i]._rewVal;
		_total_reward += _rewards_etd[i];
	}
	// _demand_pct = _dataset->_bdg_rewards_ratio[_taridx];
	// _demand = _total_reward * _demand_pct;

	fstream file(graphfile);
	if (!file) {
		cerr << "ERROR: could not open file '" << graphfile << "' for reading'" << endl;
		throw(-1);
	}

	_tmp_graph.resize(_graphsize);
	for(int i=0; i<_graphsize; i++)
		_tmp_graph[i].resize(_graphsize, false);

	string::size_type sz = 0;
	int nb_edges = -1;
	for (int i = 0; i < 6; i++){
		file >> nb_edges;
		vector<string> emptystr(nb_edges);
		for(int j = 0; j < nb_edges; j++){
			file >> emptystr[j];
			// cout << emptystr[j] << ' ';
		}
		// cout << "---------------------------" << endl;
		if(i == tar){
			for(int j = 0; j < nb_edges; j++){
				auto pos = emptystr[j].find_first_of(":");
				string temp = emptystr[j].substr(0, pos);
				int u = stoi(temp.substr(sz));
				string temp_left = emptystr[j].substr(pos+1, emptystr[j].size());
				int v = stoi(temp_left.substr(sz));
				// cout << " u, v = " << u << ", " << v << '\t';
				_tmp_graph[u][v] = true;
				_tmp_graph[v][u] = true;
			}
		}	
	}

	// exit(0);
	_graph.resize(_graphsize);
	for(int i=0; i < _graphsize; i++)
		_graph[i].resize(_graphsize, INF);
	Vertex p_i, p_j;
	for(int i=0; i < _graphsize; i++){
		for(int j=i+1; j < _graphsize; j++){
			if(_tmp_graph[i][j] == true){
				p_i._x = _dataset->_target_locs[_taridx]._x + _dataset->_all_obps[_taridx][i]._rad * cos(_dataset->_all_obps[_taridx][i]._angle);
				p_i._y = _dataset->_target_locs[_taridx]._y + _dataset->_all_obps[_taridx][i]._rad * sin(_dataset->_all_obps[_taridx][i]._angle);
				p_j._x = _dataset->_target_locs[_taridx]._x + _dataset->_all_obps[_taridx][j]._rad * cos(_dataset->_all_obps[_taridx][j]._angle);
				p_j._y = _dataset->_target_locs[_taridx]._y + _dataset->_all_obps[_taridx][j]._rad * sin(_dataset->_all_obps[_taridx][j]._angle);
				double dist = eucl_distance(p_i, p_j);
				_graph[i][j] = dist;
				_graph[j][i] = dist;	
			}	
		}
	}	
	
	_domi_labels.resize(_graphsize);
	
	_L1.resize(_graphsize,make_pair(INF, 0.0));
	_L2.resize(_graphsize,make_pair(INF, 0.0));
	_L3.resize(_graphsize,make_pair(INF, 0.0));

}


void Pulse::recursive_search(int curnode, vector<int> curpath, double pathreward, double pathrisk, ProgTime& timeout, bool log_on){
	if(_source == _sink)
		cerr << "Current code does not work when source is same as sink!! Regenerate the grpah "<< endl;
	if(true){
		timeout.end_prog();
		if(timeout._elapsed_secs >= 600.0){
			return;
		}
	}
	_iterations++;
	// if(_iterations > 50)
	// 	exit(0);
	if (log_on){
		cout << '\n';
		cout << "@@" << _iterations << ". current path: [";
		for(unsigned i =0; i < curpath.size(); i++){
			cout << curpath[i] << ',';
		}
		cout << "], path reward: " << pathreward << ", path risk: " << pathrisk <<  ", consider moving to node " << curnode <<  endl;
	}
	if(curnode == _sink ){
		if(pathreward >= _demand && pathrisk + _graph[curpath.back()][_sink] < _curbest_objval){
			_curbest_objval = pathrisk + _graph[curpath.back()][_sink];
			curpath.push_back(_sink);
			_curbest_path.clear();
			_curbest_path = curpath;
			if(log_on){
				cout << "!!!!!!!GET AN IMPROVEMENT: new obj: " << _curbest_objval << '\t' << " new path: [" << endl;
				for(unsigned i=0; i<_curbest_path.size(); i++){
					cout << _curbest_path[i] << ' ';
				}
				cout << "\n"; 
			}
		}else{
			if(log_on){
				cout << "... reach the sink but does not improve  " << endl;
			}
		}
		return;
	}


	if(curnode != _source && curnode != _sink){
		int lastnode = curpath.back();
		if(pathrisk + _lbrisk_sink[curnode] < _curbest_objval){ //primal bound check
			if(log_on)
				cout << "... still a valid partial path: check its budget " << endl;
			bool flag = is_intersected(_lbriskpaths_sink[curnode], curpath);
			if(log_on){
				cout << "curpath: ";
				for(unsigned j=0; j <curpath.size(); j++){
					cout << curpath[j] << ' ';
				}
				cout << "G(path~~>sink): ";
				for(unsigned j=0; j <_lbriskpaths_sink[curnode].size(); j++){
					cout << _lbriskpaths_sink[curnode][j] << ' ';
				}
				cout << "intersected? " << flag << '\n';
			}
			double add_reward = (flag)? 0.0:_rewards_lbrpaths_sink[curnode]-_rewards_etd[curnode];
			if(pathreward + _rewards_etd[curnode] + add_reward +0.000001 >= _demand){
				if(log_on){
					cout << "... ENOUGH BUDGET IS MET ... Pick extended path ";
				}
				pair<double,vector<int>> ret;
				if(flag == true){
					/*expand curpath to a full path by linking to shortest path to sink over subgraph */
					ret = quick_wrapup(curpath, curnode);
					if(log_on){
						cout << "on **Subgraph**: ";
						for(unsigned j=0; j <ret.second.size(); j++){
							cout << ret.second[j] << ' ';
						}
						cout << '\n';
					}
				}else{
					ret = make_pair(_lbrisk_sink[curnode], _lbriskpaths_sink[curnode]);
					if(log_on){
						cout << "on **orginal Graph**: ";
						for(unsigned j=0; j <ret.second.size(); j++){
							cout << ret.second[j] << ' ';
						}
						cout << '\n';
					}
				}

				if (pathrisk + _graph[lastnode][curnode] + ret.first < _curbest_objval){ // primal bound update 
					_curbest_objval = pathrisk + _graph[lastnode][curnode] + ret.first;
					_curbest_path.clear();
					_curbest_path.insert(_curbest_path.end(), curpath.begin(), curpath.end());
					_curbest_path.insert(_curbest_path.end(), ret.second.begin(), ret.second.end());
					if(log_on){
						cout << "!!!!!!!!find a better path. Total risk: " << _curbest_objval << " new path: ";
						for(unsigned i=0; i<_curbest_path.size(); i++){
							cout << _curbest_path[i] << ' ';
						}
						cout << "curbest_obj: " << _curbest_objval <<  "\n"; 
					}							
				}else{
					if(log_on)	
						cout << "... Propogate to sink. No better than cur OBJ." << endl;
				}
			}else{
				if(log_on){
					cout << "... NOT MEET BUDGET YET. NEED TO CONTINUE... " << endl;
				}
				auto newpath = curpath;
				newpath.push_back(curnode);
	
				double newpathrisk = pathrisk + _graph[lastnode][curnode];
				double newpathreward = pathreward + _rewards_etd[curnode];
				bool is_dominated = check_dominance(curnode, newpathrisk, newpathreward);
				if(!is_dominated){
					if(log_on)
						cout <<"... pass dominance check" << endl;
					update_domilabels(curnode, newpathrisk, newpathreward);
					for(int ngb=0; ngb <_graphsize; ngb++){
						if(_graph[curnode][ngb] < 1000.0){ //  && abs(_rewards_etd[ngb]) > 0.00001
							if(std::find(newpath.begin(), newpath.end(), ngb) == newpath.end()){			
								if(log_on){
									cout << "... " << ngb << " is not visited. *propogate* Input path: [";
									for(unsigned i=0; i<newpath.size(); i++)
										cout << newpath[i] << ' ';

									cout << ']' << endl;
								}								
								this->recursive_search(ngb, newpath, newpathreward, newpathrisk, timeout, log_on);
							}else{
								if(log_on)
									cout << "... " << ngb << " is already visited" << endl;
							}							
						}
					}
				}else{
					if(log_on)
						cout <<"... pruned by dominance" << endl;
				}
			}

		}else{
			if(log_on)
				cout << "... discarded by primal bound. risk of best full path: " << pathrisk + _lbrisk_sink[curnode] << " curbest:" << _curbest_objval << endl;
		}
	}else{ 
		/* current node is the source. Propogate to every neighbor */
		vector<int> curpath = {_source};
		for(int ngb=0; ngb<_graphsize; ngb++){
			if(_graph[curnode][ngb] < 10000.0 ){ // && abs(_rewards_etd[ngb]) > 0.00001
				// cout << "... Propogate from source " << _source << " to node " << ngb << endl;
				this->recursive_search(ngb, curpath, 0.0, 0.0, timeout, log_on);
			}
		}
	}
}

bool Pulse::check_dominance(int node, double pathrisk, double pathreward){
	set<pair<double,double>>::iterator itr;
	for(itr=_domi_labels[node].begin(); itr !=_domi_labels[node].end(); ){
		// cout << node << " cur: " << pathrisk << ", " << pathreward << " with: " <<itr->first << ", " << itr->second << endl;
		if((pathrisk >= itr->first && pathreward <= itr->second)){		
			// cout << " -------unfinished---------" << endl;
			return true;
		}
		++itr;
	}
	// cout << " -------done---------" << endl;
	return false;

	// if((pathrisk >= _L1[node].first && pathreward <= _L1[node].second)||(pathrisk >= _L2[node].first && pathreward <= _L2[node].second)||(pathrisk >= _L3[node].first && pathreward <= _L3[node].second)){
	// 	return true;
	// }
	// return false;

}

void Pulse::update_domilabels(int node, double pathrisk, double pathreward){
	if(_domi_labels[node].size() == 0){
		_domi_labels[node].insert(make_pair(pathrisk, pathreward));
	}else{
		bool IN = true;
		set<pair<double,double>>::iterator itr;
		for(itr=_domi_labels[node].begin(); itr !=_domi_labels[node].end();){
			if(pathrisk <= itr->first && pathreward >= itr->second){				
				itr = _domi_labels[node].erase(itr);
			}else{
				itr++;
			}
			if(itr->first <= pathrisk && itr->second >= pathreward){				
				IN = false;
			}
		}
		if(_domi_labels[node].size() < _maxlabsize && IN){
			_domi_labels[node].insert(make_pair(pathrisk, pathreward));
		}
	}


	// if(pathrisk <= _L1[node].first){
	// 	_L1[node] = make_pair(pathrisk, pathreward );
	// }else if(pathreward <= _L2[node].second){
	// 	_L2[node] = make_pair(pathrisk, pathreward);
	// }else{
	// 	_L3[node] = make_pair(pathrisk, pathreward);
	// }

}

bool Pulse::is_intersected(vector<int>& vec1, vector<int>& vec2){
	for(unsigned i=0; i<vec1.size(); i++){
		for(unsigned j=0; j<vec2.size(); j++){
			if(vec1[i] == vec2[j])
				return true;
		}
	}
	return false;
}


void Pulse::calc_leastrisk_sink(){
	// cout << "====>>>> Algorithm iniialization: least risky path to sink " << endl;
	vector<bool> visited(_graphsize, false);
	// int nb_voidpoints = 0;
	// for(int i=0; i<_graphsize; i++){
	// 	if(abs(_rewards_etd[i])<0.0000001){
	// 		visited[i] = true;
	// 		// nb_voidpoints++;
	// 	}
	// }
	// visited[_source] = false;
	// visited[_sink] = false;
	_lbrisk_sink.resize(_graphsize, INF);
	_lbrisk_sink[_sink] = 0.0;
	// cout << "sink: " << _sink;
	vector<int> prevnode(_graphsize, -2);


	for(int i=0; i< _graphsize; i++){
		/* find the u with smallest value dist[u]*/
		int u = -1;
		double MIN = INF;
		for(int j=0; j<_graphsize; j++){
			if(visited[j]==false && _lbrisk_sink[j] < MIN){
				u = j;
				MIN = _lbrisk_sink[j];
			}
		}
		// cout << "u = " << u << endl;
		if(u==-1){
			cout << "the orginal graph is disconnected" << endl;
			exit(0);
		}
		visited[u] = true;
		for(int v=0; v<_graphsize; v++){
			if(visited[v] == false && _lbrisk_sink[u] + _graph[u][v] < _lbrisk_sink[v]){ //*** not use void edge
				_lbrisk_sink[v] = _lbrisk_sink[u] + _graph[u][v];
				prevnode[v] = u;
			}
		}
	}
	// get the optimal path by backtracking
	_lbriskpaths_sink.resize(_graphsize);
	_rewards_lbrpaths_sink.resize(_graphsize, 0.0);
	for(int i=0; i<_graphsize; i++){
		int curnode = i;
		while(true){
			if(curnode == -2)
				break;
			_rewards_lbrpaths_sink[i] += _rewards_etd[curnode];
			_lbriskpaths_sink[i].push_back(curnode);
			curnode = prevnode[curnode];		
		}
	}
	if(false){
		for(int i=0; i<_graphsize; i++){
			cout << " ... node " << i << " to sink " << _sink << ": ";
			for(unsigned j=0; j<_lbriskpaths_sink[i].size();j++){
				cout << _lbriskpaths_sink[i][j] << ',';
			}
			cout << " ---- Total rewards: " << _rewards_lbrpaths_sink[i];
			cout << '\n';
		}
	}
	_demand =  (_total_reward/3.0 - _rewards_lbrpaths_sink[_source]) * _demand_pct + _rewards_lbrpaths_sink[_source];

	if (_demand <= _rewards_lbrpaths_sink[_source]){
		cout << "Reward on min-risk path already satisfies the demand!!!" << endl;
		// exit(0);
	}
}

pair<double, vector<int>> Pulse::quick_wrapup(const vector<int>& visited, int startnode){
	vector<bool> visited_etd(_graphsize, false);
	for(unsigned i=0; i<visited.size(); i++)
		visited_etd[visited[i]] = true;
	// for(unsigned i = 0; i < visited.size(); i++)
	// 	cout << visited[i] << ' ';
	// cout << '\n';
	vector<int> prevnode(_graphsize, -2);
	vector<double> dist(_graphsize, INF);
	dist[startnode] = 0.0;
	for(int i=0; i< _graphsize-(int)visited.size(); i++){
		/* find the u with smallest value dist[u]*/
		int u = -1;
		double MIN = INF;
		for(int j=0; j<_graphsize; j++){
			if(visited_etd[j]==false && dist[j] < MIN){
				u = j;
				MIN = dist[j];
			}
		}
		if(u==-1 && visited_etd[_sink] == false){
			// cerr << "the orginal graph is disconnected" << endl;
			break;
		}
		visited_etd[u] = true;
		if(u == _sink){
			break;
		}
		for(int v=0; v<_graphsize; v++){
			if(visited_etd[v] == false && dist[u] + _graph[u][v] < dist[v]){
				dist[v] = dist[u] + _graph[u][v];
				prevnode[v] = u;
			}
		}


	}

	// cout << "prevnode: ";
	// for(int i=0; i<_graphsize; i++)
	// 	cout << prevnode[i] << ' ';
	// cout << '\n';
	// cout << "dist: ";
	// for(int i=0; i<_graphsize; i++)
	// 	cout << dist[i] << ' ';
	// cout << '\n';
	vector<int> optpath;
	if(visited_etd[_sink]){
		// get the optimal path by backtracking
		int curnode = _sink;
		while(curnode!= -2){
			optpath.insert(optpath.begin(),curnode);		
			curnode = prevnode[curnode];
		}
	}
	// cout << "path~~>sink: ";
	// for(unsigned i=0; i<optpath.size(); i++)
	// 	cout << optpath[i] << ' ';
	// cout << '\n';
	return make_pair(dist[_sink], optpath);
}



void Pulse::initialize_generalgraph(vector<vector<double>>& adjmatx, vector<double>& rewards, int size, double budget_pct){
	this->_graphsize = size;
	this->_source = 0;
	this->_sink = size -1;

	_rewards_etd = rewards;
	_total_reward = 0.0;
	for(int i=0; i< size; i++){
		// cout << _rewards_etd[i] << " ";
		_total_reward += _rewards_etd[i];
	}
	// cout << "total " << _total_reward << endl;
	_demand_pct = budget_pct;
	_demand = _total_reward * _demand_pct;
	_graph = adjmatx;

	for(int i=0; i<size; i++){
		for(int j=0; j<size; j++){
			if (_graph[i][j] == 0){
				_graph[i][j] = INF;
			}
			// cout << _graph[i][j] << " ";
		}
		// cout << '\n';
	}
	// _domi_labels.resize(_graphsize, make_pair(INF, 0));
	
}


double Pulse::eucl_distance(const Vertex& p1, const Vertex& p2){
	return sqrt(pow(p1._x -p2._x,2) + pow(p1._y - p2._y,2));
}



double Pulse::calc_path_rewards(vector<int> & path){
	double total_reward = 0.0;
	for(unsigned i=0; i < path.size(); i++){
		total_reward += _rewards_etd[path[i]];
	}
	return total_reward;
}


void Pulse::print_opt_sol(){
	cout << " =====>>> total reward: " << _total_reward << ", Demand: " << _demand << ", Reward on min-risk path: " << _rewards_lbrpaths_sink[_source] << endl;
	cout << " =====>>> optimal obj value: " << _curbest_objval << ", collected reward: " << calc_path_rewards(_curbest_path) << endl; 
	cout << "          optimal path [";
	double verified_obj = 0;
	for(unsigned i =0; i < _curbest_path.size(); i++){
		cout << _curbest_path[i] << ' ';
		if(i < _curbest_path.size()-1)
			verified_obj += _graph[_curbest_path[i]][_curbest_path[i+1]];
	}
	cout << "]\n";
}


// void Pulse::write_opt_sol(){
// 	ofstream file;
// 	string cur_dir = "/projects/academic/josewalt/caigao/RRARP_LinKern/ret/pulse_out/";
//     file.open(cur_dir + name_only + "_algo_" + to_string(algo_idx) + ".FisTest_" + to_string(fischetti_on));
// 	file << "---> instance_name: " << instance << '\n';
// 	file.close();
// }