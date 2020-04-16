#include "PulseAlgo.h"
#define _USE_MATH_DEFINES
#include<bits/stdc++.h>


void Pulse::initialize(DataHandler& dataset, int tar, int nb_dstzn, int lmt_obps, string graphfile, int source, int sink){
	this->_dataset = &dataset;
	this->_taridx = tar;
	this->_graphsize = lmt_obps + 16;
	this->_source = source;
	this->_sink = sink;

	_rewards_etd.resize(_graphsize, 0);
	_total_reward = 0.0;
	for(int i=0; i< _graphsize; i++){
		_rewards_etd[i] = _dataset->_all_obps[_taridx][i]._rewVal;
		_total_reward += _rewards_etd[i];
		// if(_rewards_etd[i] > 0)
		// 	cout << _rewards_etd[i] << " ";
	}

	// cout << '\n';
	// cout << _total_reward << endl;

	_budget_pct = _dataset->_bdg_rewards_ratio[_taridx];
	_budget_pct = 0.1;

	_budget = _total_reward * _budget_pct;
	// cout << "total reward: " << _total_reward << '\t' ;
	// cout << "_budget: " << _budget << endl;

	fstream file(graphfile);
	if (!file) {
		cerr << "ERROR: could not open file '" << graphfile << "' for reading'" << endl;
		throw(-1);
	}

	_tmp_graph.resize(_graphsize);
	for(int i=0; i<_graphsize; i++)
		_tmp_graph[i].resize(_graphsize, false);

	string::size_type sz;
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
	
	_domi_labels.resize(_graphsize, make_pair(INF, 0));

	// for(int i=0; i < _graphsize ; i++){
	// 	for(int j=0; j < _graphsize; j++){
	// 		cout << _tmp_graph[i][j] << '\t';
	// 	}	
	// 	cout << '\n';
	// }

}














void Pulse::initialize(DataHandler& dataset, int tar, int nb_dstzn, int lmt_obps){
	this->_dataset = &dataset;
	this->_taridx = tar;
	this->_graphsize = lmt_obps + 2;
	this->_source = 0;
	this->_sink = lmt_obps + 1;

	_rewards_etd.resize(_graphsize, 0);
	_total_reward = 0.0;
	for(int i=1; i<= lmt_obps; i++){
		_rewards_etd[i] = _dataset->_all_obps[_taridx][i-1]._rewVal;
		_total_reward += _rewards_etd[i];
		// cout << _rewards_etd[i] << " ";
	}

	// cout << '\n';
	// cout << _total_reward << endl;

	_budget_pct = _dataset->_bdg_rewards_ratio[_taridx];
	_budget_pct = 0.5;

	_budget = _total_reward * _budget_pct;
	// cout << "total reward: " << _total_reward << '\t' ;
	// cout << "_budget: " << _budget << endl;

	_graph.resize(_graphsize);
	for(int i=0; i < _graphsize; i++)
		_graph[i].resize(_graphsize, INF);
	Vertex p_i, p_j;
	for(int i=1; i <= lmt_obps; i++){
		for(int j=i+1; j <= lmt_obps; j++){
			p_i._x = _dataset->_target_locs[_taridx]._x + _dataset->_all_obps[_taridx][i-1]._rad * cos(_dataset->_all_obps[_taridx][i-1]._angle);
			p_i._y = _dataset->_target_locs[_taridx]._y + _dataset->_all_obps[_taridx][i-1]._rad * sin(_dataset->_all_obps[_taridx][i-1]._angle);
			p_j._x = _dataset->_target_locs[_taridx]._x + _dataset->_all_obps[_taridx][j-1]._rad * cos(_dataset->_all_obps[_taridx][j-1]._angle);
			p_j._y = _dataset->_target_locs[_taridx]._y + _dataset->_all_obps[_taridx][j-1]._rad * sin(_dataset->_all_obps[_taridx][j-1]._angle);
			double dist = eucl_distance(p_i, p_j);
			if(dist < 0.5){
				_graph[i][j] = dist;
				_graph[j][i] = dist;	
			}	
		}	
	}

	this->_nb_tnps = nb_dstzn;
	_tnps_tar.resize(_nb_tnps);
	double subarc_angle = 2.0 * M_PI /(double)_nb_tnps;
	for(int i =0; i < _nb_tnps; i++){
		_tnps_tar[i]._x = _dataset->_radii[_taridx] * cos(i * subarc_angle) + _dataset->_target_locs[_taridx]._x;
		_tnps_tar[i]._y = _dataset->_radii[_taridx] * sin(i * subarc_angle) + _dataset->_target_locs[_taridx]._y;
		// _tnps_tar[i].print();
	}

	_domi_labels.resize(_graphsize, make_pair(INF, 0));

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
	_budget_pct = budget_pct;
	_budget = _total_reward * _budget_pct;
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
	_domi_labels.resize(_graphsize, make_pair(INF, 0));
	
}

void Pulse::update_graph(Vertex& entry, Vertex& exit){
		// source and sink
	// cout << "entry and exit :" << endl;
	// entry.print();
	// exit.print();
	Vertex p;
	for(int i=1; i <= _graphsize-2; i++){
		p._x = _dataset->_target_locs[_taridx]._x + _dataset->_all_obps[_taridx][i-1]._rad * cos(_dataset->_all_obps[_taridx][i-1]._angle);
		p._y = _dataset->_target_locs[_taridx]._y + _dataset->_all_obps[_taridx][i-1]._rad * sin(_dataset->_all_obps[_taridx][i-1]._angle);
		// p.print();
		double dist = eucl_distance(entry, p);
		if(dist < 0.3){
			_graph[_source][i] = dist;
			_graph[i][_source] = dist;
		}else{
			_graph[_source][i] = INF;
			_graph[i][_source] = INF;
		}
		dist = eucl_distance(exit, p);
		if(dist < 0.3){
			_graph[_sink][i] = dist;
			_graph[i][_sink] = dist;
		}else{
			_graph[_sink][i] = INF;
			_graph[i][_sink] = INF;
		}
	}


	// for(int i=0; i < _graphsize ; i++){
	// 	for(int j=0; j < _graphsize; j++){
	// 		cout << _graph[i][j] << '\t';
	// 	}	
	// 	cout << '\n';
	// }
}

double Pulse::eucl_distance(const Vertex& p1, const Vertex& p2){
	return sqrt(pow(p1._x -p2._x,2) + pow(p1._y - p2._y,2));
}

void Pulse::recursive_search(int curnode, vector<int> curpath, double pathreward, double pathrisk, bool log_on){
	if(_source == _sink)
		cerr << "Current code does not work when source is same as sink!! Regenerate the grpah "<< endl;
	_iterations++;
	if (log_on){
		cout << '\n';
		cout << "@@" << _iterations << ". current path: [";
		for(unsigned i =0; i < curpath.size(); i++){
			cout << curpath[i] << ',';
		}
		cout << "], path reward: " << pathreward << ", path risk: " << pathrisk <<  ", consider moving to node " << curnode <<  endl;
	}
	if(curnode == _sink ){
		if(pathreward >= _budget && pathrisk + _graph[curpath.back()][_sink] < _curbest_objval){
			_curbest_objval = pathrisk + _graph[curpath.back()][_sink];
			curpath.push_back(_sink);
			_curbest_path.clear();
			_curbest_path = curpath;
			// cout << "obj* (sink), " << _curbest_objval << " @iterations " << _iterations << endl;

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


	if(curnode != _source){
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
			if(pathreward + _rewards_etd[curnode] + add_reward >= _budget){// partial path expansion 
				if(log_on){
					cout << "... ENOUGH BUDGET IS MET ..." ;
					if(flag){
						cout << " expanded by (subgraph(v~~>b)); curpath length: " << curpath.size() << endl; 
					}else{
						cout << " expanded by (orginal graph (v~~>b)); curpath length: " << curpath.size()  << endl;
					}

				}
				pair<double,vector<int>> ret;
				if(flag == true){
					/*expand curpath to a full path by linking to shortest path to sink over subgraph */
					ret = quick_wrapup(curpath, curnode);
					if(log_on){
						cout << "path_to_sink(Subgraph): ";
						for(unsigned j=0; j <ret.second.size(); j++){
							cout << ret.second[j] << ' ';
						}
						cout << '\n';
					}
				}else{
					ret = make_pair(_lbrisk_sink[curnode], _lbriskpaths_sink[curnode]);
					if(log_on){
						cout << "path_to_sink(orginal Graph): ";
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
				}
			}else{
				if(log_on){
					cout << "... NOT MEET BUDGET YET. NEED TO CONTINUE... " << endl;
				}
				auto newpath = curpath;
				newpath.push_back(curnode);
	
				double newpathrisk = pathrisk + _graph[lastnode][curnode];
				double newpathreward = pathreward + _rewards_etd[curnode];
				if(!check_dominance(curnode, newpathrisk, newpathreward)){
					if(log_on)
						cout <<"... pass dominance check" << endl;
					update_domilabels(curnode, newpathrisk, newpathreward);
					for(int ngb=0; ngb <_graphsize; ngb++){
						if(_graph[curnode][ngb] < 1000.0 && abs(_rewards_etd[ngb]) > 0.00001){
							if(std::find(newpath.begin(), newpath.end(), ngb) == newpath.end()){			
								if(log_on){
									cout << "... " << ngb << " is not visited. Move to it. input path: [";
									for(unsigned i=0; i<newpath.size(); i++)
										cout << newpath[i] << ' ';

									cout << ']' << endl;
								}								
								this->recursive_search(ngb, newpath, newpathreward, newpathrisk, log_on);
							}else{
								if(log_on)
									cout << "... " << ngb << " is already visited" << endl;
							}							
						}
					}
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
			if(_graph[curnode][ngb] < 10000.0 && abs(_rewards_etd[ngb]) > 0.00001){
				cout << "... Propogate from source " << _source << " to node " << ngb << endl;
				this->recursive_search(ngb, curpath, 0.0, 0.0, log_on);
			}
		}
	}
}

bool Pulse::check_dominance(int node, double pathrisk, double pathreward){
	if( (pathrisk >= _domi_labels[node].first && pathreward < _domi_labels[node].second) || (pathrisk > _domi_labels[node].first && pathreward <= _domi_labels[node].second)){
		return true;
	}else{
		return false;
	}

}

void Pulse::update_domilabels(int node, double pathrisk, double pathreward){
	if(pathrisk <= _domi_labels[node].first && pathreward >= _domi_labels[node].second){
		_domi_labels[node].first = pathrisk;
		_domi_labels[node].second = pathreward;
	}
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
	cout << "====>>>> Algorithm iniialization: least risky path to sink " << endl;
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
			cerr << "the orginal graph is disconnected" << endl;
		}
		visited[u] = true;
		for(int v=0; v<_graphsize; v++){
			if(visited[v] == false && _lbrisk_sink[u] + _graph[u][v] < _lbrisk_sink[v] && !(u<16 && v<16)){ //*** not use void edge
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
	for(int i=0; i<_graphsize; i++){
		cout << " ... node " << i << " to sink " << _sink << ": ";
		for(unsigned j=0; j<_lbriskpaths_sink[i].size();j++){
			cout << _lbriskpaths_sink[i][j] << ',';
		}
		cout << " ---- Total rewards: " << _rewards_lbrpaths_sink[i];
		cout << '\n';
	}
	// // exit(1);

	// cout << "least risk path to sink: " ;
	// for(int i=0; i<_graphsize; i++)
	// 	cout << _lbrisk_sink[i] << " ";
	// cout << endl;
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

void Pulse::print_opt_sol(){
	cout << " =====>>> total reward: " << _total_reward << ", required budget: " << _budget << endl;
	cout << " =====>>> optimal obj value: " << _curbest_objval << "  optimal path: ";
	for(unsigned i =0; i < _curbest_path.size(); i++)
		cout << _curbest_path[i] << ' ';
	cout << '\n';
}
